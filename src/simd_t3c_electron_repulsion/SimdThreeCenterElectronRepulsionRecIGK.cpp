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


#include "SimdThreeCenterElectronRepulsionRecIGK.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIF.hpp"
#include "SimdTransferIG.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKF.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLD.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransferMP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_igk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_igk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 214487, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 214487, 165887, 10740, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto mu = a_exps[i] * b_exps[j] / p;

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fa = -b_exps[j] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pa(buffer, coordinates, 0, nmax, fa);

                simdfunc::compute_pc(buffer, coordinates, c_coordinates, 3, n, nmax, fc);

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 6, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                                        16, 17}, ncols, fj, mu, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 24, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 27, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 30, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 33, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 36, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 39, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 42, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 45, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 48, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 51, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 54, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 57, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 60, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 63, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 66, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 69, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 72, 0, 3, 7, 8,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 78, 0, 3, 8, 9,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 84, 0, 3, 9, 10,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 90, 0, 3, 10, 11,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 96, 0, 3, 11, 12,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 102, 0, 3, 12, 13,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 108, 0, 3, 13, 14,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 114, 0, 3, 14, 15,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 120, 0, 3, 15, 16,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 126, 0, 3, 16, 17,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 132, 0, 3, 17, 18,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 138, 0, 3, 18, 19,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 144, 0, 3, 19, 20,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 150, 0, 3, 20, 21,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 156, 0, 3, 21, 22,
                                                                       66, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 162, 0, 3, 24, 27,
                                                                       72, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 172, 0, 3, 27, 30,
                                                                       78, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 182, 0, 3, 30, 33,
                                                                       84, 90, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 192, 0, 3, 33, 36,
                                                                       90, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 202, 0, 3, 36, 39,
                                                                       96, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 212, 0, 3, 39, 42,
                                                                       102, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 222, 0, 3, 42, 45,
                                                                       108, 114, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 232, 0, 3, 45, 48,
                                                                       114, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 242, 0, 3, 48, 51,
                                                                       120, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 252, 0, 3, 51, 54,
                                                                       126, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 262, 0, 3, 54, 57,
                                                                       132, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 272, 0, 3, 57, 60,
                                                                       138, 144, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 282, 0, 3, 60, 63,
                                                                       144, 150, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 292, 0, 3, 63, 66,
                                                                       150, 156, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 302, 0, 3, 72, 78,
                                                                       162, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 317, 0, 3, 78, 84,
                                                                       172, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 332, 0, 3, 84, 90,
                                                                       182, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 347, 0, 3, 90, 96,
                                                                       192, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 362, 0, 3, 96,
                                                                       102, 202, 212, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 377, 0, 3, 102,
                                                                       108, 212, 222, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 392, 0, 3, 108,
                                                                       114, 222, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 407, 0, 3, 114,
                                                                       120, 232, 242, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 422, 0, 3, 120,
                                                                       126, 242, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 437, 0, 3, 126,
                                                                       132, 252, 262, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 452, 0, 3, 132,
                                                                       138, 262, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 467, 0, 3, 138,
                                                                       144, 272, 282, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 482, 0, 3, 144,
                                                                       150, 282, 292, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 497, 0, 3, 162,
                                                                       172, 302, 317, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 518, 0, 3, 172,
                                                                       182, 317, 332, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 539, 0, 3, 182,
                                                                       192, 332, 347, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 560, 0, 3, 192,
                                                                       202, 347, 362, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 581, 0, 3, 202,
                                                                       212, 362, 377, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 602, 0, 3, 212,
                                                                       222, 377, 392, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 623, 0, 3, 222,
                                                                       232, 392, 407, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 644, 0, 3, 232,
                                                                       242, 407, 422, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 665, 0, 3, 242,
                                                                       252, 422, 437, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 686, 0, 3, 252,
                                                                       262, 437, 452, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 707, 0, 3, 262,
                                                                       272, 452, 467, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 728, 0, 3, 272,
                                                                       282, 467, 482, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 749, 0, 3, 302,
                                                                       317, 497, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 777, 0, 3, 317,
                                                                       332, 518, 539, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 805, 0, 3, 332,
                                                                       347, 539, 560, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 833, 0, 3, 347,
                                                                       362, 560, 581, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 861, 0, 3, 362,
                                                                       377, 581, 602, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 889, 0, 3, 377,
                                                                       392, 602, 623, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 917, 0, 3, 392,
                                                                       407, 623, 644, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 945, 0, 3, 407,
                                                                       422, 644, 665, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 973, 0, 3, 422,
                                                                       437, 665, 686, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1001, 0, 3, 437,
                                                                       452, 686, 707, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1029, 0, 3, 452,
                                                                       467, 707, 728, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1057, 0, 3, 497,
                                                                       518, 749, 777, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1093, 0, 3, 518,
                                                                       539, 777, 805, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1129, 0, 3, 539,
                                                                       560, 805, 833, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1165, 0, 3, 560,
                                                                       581, 833, 861, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1201, 0, 3, 581,
                                                                       602, 861, 889, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1237, 0, 3, 602,
                                                                       623, 889, 917, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1273, 0, 3, 623,
                                                                       644, 917, 945, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1309, 0, 3, 644,
                                                                       665, 945, 973, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1345, 0, 3, 665,
                                                                       686, 973, 1001, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1381, 0, 3, 686,
                                                                       707, 1001, 1029, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1417, 0, 3, 749,
                                                                       777, 1057, 1093, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1462, 0, 3, 777,
                                                                       805, 1093, 1129, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1507, 0, 3, 805,
                                                                       833, 1129, 1165, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1552, 0, 3, 833,
                                                                       861, 1165, 1201, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1597, 0, 3, 861,
                                                                       889, 1201, 1237, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1642, 0, 3, 889,
                                                                       917, 1237, 1273, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1687, 0, 3, 917,
                                                                       945, 1273, 1309, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1732, 0, 3, 945,
                                                                       973, 1309, 1345, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1777, 0, 3, 973,
                                                                       1001, 1345, 1381, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1822, 0, 3, 1057,
                                                                       1093, 1417, 1462, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1877, 0, 3, 1093,
                                                                       1129, 1462, 1507, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1932, 0, 3, 1129,
                                                                       1165, 1507, 1552, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1987, 0, 3, 1165,
                                                                       1201, 1552, 1597, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2042, 0, 3, 1201,
                                                                       1237, 1597, 1642, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2097, 0, 3, 1237,
                                                                       1273, 1642, 1687, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2152, 0, 3, 1273,
                                                                       1309, 1687, 1732, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2207, 0, 3, 1309,
                                                                       1345, 1732, 1777, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2262, 0, 3, 1417,
                                                                       1462, 1822, 1877, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2328, 0, 3, 1462,
                                                                       1507, 1877, 1932, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2394, 0, 3, 1507,
                                                                       1552, 1932, 1987, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2460, 0, 3, 1552,
                                                                       1597, 1987, 2042, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2526, 0, 3, 1597,
                                                                       1642, 2042, 2097, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2592, 0, 3, 1642,
                                                                       1687, 2097, 2152, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2658, 0, 3, 1687,
                                                                       1732, 2152, 2207, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2724, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2727, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2730, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2733, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2736, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2739, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2742, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2745, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2748, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2751, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2754, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2757, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2760, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2763, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2766, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2769, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2772, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2775, 3, 7, 24,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2784, 3, 8, 27,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2793, 3, 9, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2802, 3, 10, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2811, 3, 11, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2820, 3, 12, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2829, 3, 13, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2838, 3, 14, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2847, 3, 15, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2856, 3, 16, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2865, 3, 17, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2874, 3, 18, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2883, 3, 19, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2892, 3, 20, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2901, 3, 21, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2910, 3, 22, 69,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2919, 3, 24, 72,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2937, 3, 27, 78,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2955, 3, 30, 84,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2973, 3, 33, 90,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2991, 3, 36, 96,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3009, 3, 39, 102,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3027, 3, 42, 108,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3045, 3, 45, 114,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3063, 3, 48, 120,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3081, 3, 51, 126,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3099, 3, 54, 132,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3117, 3, 57, 138,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3135, 3, 60, 144,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3153, 3, 63, 150,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3171, 3, 66, 156,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3189, 3, 72, 162,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3219, 3, 78, 172,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3249, 3, 84, 182,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3279, 3, 90, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3309, 3, 96, 202,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3339, 3, 102, 212,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3369, 3, 108, 222,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3399, 3, 114, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3429, 3, 120, 242,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3459, 3, 126, 252,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3489, 3, 132, 262,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3519, 3, 138, 272,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3549, 3, 144, 282,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3579, 3, 150, 292,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3609, 3, 162, 302,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3654, 3, 172, 317,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3699, 3, 182, 332,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3744, 3, 192, 347,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3789, 3, 202, 362,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3834, 3, 212, 377,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3879, 3, 222, 392,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3924, 3, 232, 407,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3969, 3, 242, 422,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4014, 3, 252, 437,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4059, 3, 262, 452,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4104, 3, 272, 467,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4149, 3, 282, 482,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4194, 3, 302, 497,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4257, 3, 317, 518,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4320, 3, 332, 539,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4383, 3, 347, 560,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4446, 3, 362, 581,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4509, 3, 377, 602,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4572, 3, 392, 623,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4635, 3, 407, 644,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4698, 3, 422, 665,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4761, 3, 437, 686,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4824, 3, 452, 707,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4887, 3, 467, 728,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4950, 3, 497, 749,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5034, 3, 518, 777,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5118, 3, 539, 805,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5202, 3, 560, 833,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5286, 3, 581, 861,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5370, 3, 602, 889,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5454, 3, 623, 917,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5538, 3, 644, 945,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5622, 3, 665, 973,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5706, 3, 686,
                                                                       1001, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5790, 3, 707,
                                                                       1029, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5874, 3, 749,
                                                                       1057, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5982, 3, 777,
                                                                       1093, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6090, 3, 805,
                                                                       1129, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6198, 3, 833,
                                                                       1165, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6306, 3, 861,
                                                                       1201, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6414, 3, 889,
                                                                       1237, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6522, 3, 917,
                                                                       1273, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6630, 3, 945,
                                                                       1309, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6738, 3, 973,
                                                                       1345, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6846, 3, 1001,
                                                                       1381, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6954, 3, 1057,
                                                                       1417, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7089, 3, 1093,
                                                                       1462, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7224, 3, 1129,
                                                                       1507, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7359, 3, 1165,
                                                                       1552, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7494, 3, 1201,
                                                                       1597, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7629, 3, 1237,
                                                                       1642, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7764, 3, 1273,
                                                                       1687, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7899, 3, 1309,
                                                                       1732, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8034, 3, 1345,
                                                                       1777, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8169, 3, 1417,
                                                                       1822, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8334, 3, 1462,
                                                                       1877, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8499, 3, 1507,
                                                                       1932, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8664, 3, 1552,
                                                                       1987, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8829, 3, 1597,
                                                                       2042, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8994, 3, 1642,
                                                                       2097, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9159, 3, 1687,
                                                                       2152, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9324, 3, 1732,
                                                                       2207, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9489, 3, 1822,
                                                                       2262, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9687, 3, 1877,
                                                                       2328, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9885, 3, 1932,
                                                                       2394, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10083, 3, 1987,
                                                                       2460, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10281, 3, 2042,
                                                                       2526, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10479, 3, 2097,
                                                                       2592, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10677, 3, 2152,
                                                                       2658, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10875, 3, 7, 8,
                                                                       2730, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10881, 3, 8, 9,
                                                                       2733, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10887, 3, 9, 10,
                                                                       2736, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10893, 3, 10, 11,
                                                                       2739, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10899, 3, 11, 12,
                                                                       2742, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10905, 3, 12, 13,
                                                                       2745, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10911, 3, 13, 14,
                                                                       2748, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10917, 3, 14, 15,
                                                                       2751, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10923, 3, 15, 16,
                                                                       2754, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10929, 3, 16, 17,
                                                                       2757, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10935, 3, 17, 18,
                                                                       2760, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10941, 3, 18, 19,
                                                                       2763, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10947, 3, 19, 20,
                                                                       2766, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10953, 3, 20, 21,
                                                                       2769, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10959, 3, 21, 22,
                                                                       2772, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10965, 0, 3,
                                                                       10875, 2730, 10881, 2793,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10983, 0, 3,
                                                                       10881, 2733, 10887, 2802,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11001, 0, 3,
                                                                       10887, 2736, 10893, 2811,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11019, 0, 3,
                                                                       10893, 2739, 10899, 2820,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11037, 0, 3,
                                                                       10899, 2742, 10905, 2829,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11055, 0, 3,
                                                                       10905, 2745, 10911, 2838,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11073, 0, 3,
                                                                       10911, 2748, 10917, 2847,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11091, 0, 3,
                                                                       10917, 2751, 10923, 2856,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11109, 0, 3,
                                                                       10923, 2754, 10929, 2865,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11127, 0, 3,
                                                                       10929, 2757, 10935, 2874,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11145, 0, 3,
                                                                       10935, 2760, 10941, 2883,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11163, 0, 3,
                                                                       10941, 2763, 10947, 2892,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11181, 0, 3,
                                                                       10947, 2766, 10953, 2901,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11199, 0, 3,
                                                                       10953, 2769, 10959, 2910,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11217, 0, 3,
                                                                       10965, 2793, 10983, 72,
                                                                       78, 2955, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11253, 0, 3,
                                                                       10983, 2802, 11001, 78,
                                                                       84, 2973, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11289, 0, 3,
                                                                       11001, 2811, 11019, 84,
                                                                       90, 2991, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11325, 0, 3,
                                                                       11019, 2820, 11037, 90,
                                                                       96, 3009, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11361, 0, 3,
                                                                       11037, 2829, 11055, 96,
                                                                       102, 3027, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11397, 0, 3,
                                                                       11055, 2838, 11073, 102,
                                                                       108, 3045, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11433, 0, 3,
                                                                       11073, 2847, 11091, 108,
                                                                       114, 3063, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11469, 0, 3,
                                                                       11091, 2856, 11109, 114,
                                                                       120, 3081, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11505, 0, 3,
                                                                       11109, 2865, 11127, 120,
                                                                       126, 3099, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11541, 0, 3,
                                                                       11127, 2874, 11145, 126,
                                                                       132, 3117, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11577, 0, 3,
                                                                       11145, 2883, 11163, 132,
                                                                       138, 3135, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11613, 0, 3,
                                                                       11163, 2892, 11181, 138,
                                                                       144, 3153, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11649, 0, 3,
                                                                       11181, 2901, 11199, 144,
                                                                       150, 3171, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11685, 0, 3,
                                                                       11217, 2955, 11253, 162,
                                                                       172, 3249, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11745, 0, 3,
                                                                       11253, 2973, 11289, 172,
                                                                       182, 3279, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11805, 0, 3,
                                                                       11289, 2991, 11325, 182,
                                                                       192, 3309, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11865, 0, 3,
                                                                       11325, 3009, 11361, 192,
                                                                       202, 3339, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11925, 0, 3,
                                                                       11361, 3027, 11397, 202,
                                                                       212, 3369, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11985, 0, 3,
                                                                       11397, 3045, 11433, 212,
                                                                       222, 3399, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12045, 0, 3,
                                                                       11433, 3063, 11469, 222,
                                                                       232, 3429, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12105, 0, 3,
                                                                       11469, 3081, 11505, 232,
                                                                       242, 3459, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12165, 0, 3,
                                                                       11505, 3099, 11541, 242,
                                                                       252, 3489, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12225, 0, 3,
                                                                       11541, 3117, 11577, 252,
                                                                       262, 3519, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12285, 0, 3,
                                                                       11577, 3135, 11613, 262,
                                                                       272, 3549, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12345, 0, 3,
                                                                       11613, 3153, 11649, 272,
                                                                       282, 3579, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12405, 0, 3,
                                                                       11685, 3249, 11745, 302,
                                                                       317, 3699, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12495, 0, 3,
                                                                       11745, 3279, 11805, 317,
                                                                       332, 3744, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12585, 0, 3,
                                                                       11805, 3309, 11865, 332,
                                                                       347, 3789, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12675, 0, 3,
                                                                       11865, 3339, 11925, 347,
                                                                       362, 3834, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12765, 0, 3,
                                                                       11925, 3369, 11985, 362,
                                                                       377, 3879, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12855, 0, 3,
                                                                       11985, 3399, 12045, 377,
                                                                       392, 3924, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12945, 0, 3,
                                                                       12045, 3429, 12105, 392,
                                                                       407, 3969, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13035, 0, 3,
                                                                       12105, 3459, 12165, 407,
                                                                       422, 4014, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13125, 0, 3,
                                                                       12165, 3489, 12225, 422,
                                                                       437, 4059, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13215, 0, 3,
                                                                       12225, 3519, 12285, 437,
                                                                       452, 4104, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13305, 0, 3,
                                                                       12285, 3549, 12345, 452,
                                                                       467, 4149, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13395, 0, 3,
                                                                       12405, 3699, 12495, 497,
                                                                       518, 4320, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13521, 0, 3,
                                                                       12495, 3744, 12585, 518,
                                                                       539, 4383, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13647, 0, 3,
                                                                       12585, 3789, 12675, 539,
                                                                       560, 4446, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13773, 0, 3,
                                                                       12675, 3834, 12765, 560,
                                                                       581, 4509, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13899, 0, 3,
                                                                       12765, 3879, 12855, 581,
                                                                       602, 4572, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14025, 0, 3,
                                                                       12855, 3924, 12945, 602,
                                                                       623, 4635, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14151, 0, 3,
                                                                       12945, 3969, 13035, 623,
                                                                       644, 4698, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14277, 0, 3,
                                                                       13035, 4014, 13125, 644,
                                                                       665, 4761, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14403, 0, 3,
                                                                       13125, 4059, 13215, 665,
                                                                       686, 4824, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14529, 0, 3,
                                                                       13215, 4104, 13305, 686,
                                                                       707, 4887, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14655, 0, 3,
                                                                       13395, 4320, 13521, 749,
                                                                       777, 5118, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14823, 0, 3,
                                                                       13521, 4383, 13647, 777,
                                                                       805, 5202, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14991, 0, 3,
                                                                       13647, 4446, 13773, 805,
                                                                       833, 5286, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15159, 0, 3,
                                                                       13773, 4509, 13899, 833,
                                                                       861, 5370, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15327, 0, 3,
                                                                       13899, 4572, 14025, 861,
                                                                       889, 5454, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15495, 0, 3,
                                                                       14025, 4635, 14151, 889,
                                                                       917, 5538, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15663, 0, 3,
                                                                       14151, 4698, 14277, 917,
                                                                       945, 5622, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15831, 0, 3,
                                                                       14277, 4761, 14403, 945,
                                                                       973, 5706, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15999, 0, 3,
                                                                       14403, 4824, 14529, 973,
                                                                       1001, 5790, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16167, 0, 3,
                                                                       14655, 5118, 14823, 1057,
                                                                       1093, 6090, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16383, 0, 3,
                                                                       14823, 5202, 14991, 1093,
                                                                       1129, 6198, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16599, 0, 3,
                                                                       14991, 5286, 15159, 1129,
                                                                       1165, 6306, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16815, 0, 3,
                                                                       15159, 5370, 15327, 1165,
                                                                       1201, 6414, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17031, 0, 3,
                                                                       15327, 5454, 15495, 1201,
                                                                       1237, 6522, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17247, 0, 3,
                                                                       15495, 5538, 15663, 1237,
                                                                       1273, 6630, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17463, 0, 3,
                                                                       15663, 5622, 15831, 1273,
                                                                       1309, 6738, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17679, 0, 3,
                                                                       15831, 5706, 15999, 1309,
                                                                       1345, 6846, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17895, 0, 3,
                                                                       16167, 6090, 16383, 1417,
                                                                       1462, 7224, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18165, 0, 3,
                                                                       16383, 6198, 16599, 1462,
                                                                       1507, 7359, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18435, 0, 3,
                                                                       16599, 6306, 16815, 1507,
                                                                       1552, 7494, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18705, 0, 3,
                                                                       16815, 6414, 17031, 1552,
                                                                       1597, 7629, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18975, 0, 3,
                                                                       17031, 6522, 17247, 1597,
                                                                       1642, 7764, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 19245, 0, 3,
                                                                       17247, 6630, 17463, 1642,
                                                                       1687, 7899, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 19515, 0, 3,
                                                                       17463, 6738, 17679, 1687,
                                                                       1732, 8034, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 19785, 0, 3,
                                                                       17895, 7224, 18165, 1822,
                                                                       1877, 8499, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 20115, 0, 3,
                                                                       18165, 7359, 18435, 1877,
                                                                       1932, 8664, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 20445, 0, 3,
                                                                       18435, 7494, 18705, 1932,
                                                                       1987, 8829, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 20775, 0, 3,
                                                                       18705, 7629, 18975, 1987,
                                                                       2042, 8994, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 21105, 0, 3,
                                                                       18975, 7764, 19245, 2042,
                                                                       2097, 9159, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 21435, 0, 3,
                                                                       19245, 7899, 19515, 2097,
                                                                       2152, 9324, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 21765, 0, 3,
                                                                       19785, 8499, 20115, 2262,
                                                                       2328, 9885, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 22161, 0, 3,
                                                                       20115, 8664, 20445, 2328,
                                                                       2394, 10083, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 22557, 0, 3,
                                                                       20445, 8829, 20775, 2394,
                                                                       2460, 10281, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 22953, 0, 3,
                                                                       20775, 8994, 21105, 2460,
                                                                       2526, 10479, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 23349, 0, 3,
                                                                       21105, 9159, 21435, 2526,
                                                                       2592, 10677, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23745, 3, 2724,
                                                                       2727, 10875, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23755, 3, 2727,
                                                                       2730, 10881, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23765, 3, 2730,
                                                                       2733, 10887, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23775, 3, 2733,
                                                                       2736, 10893, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23785, 3, 2736,
                                                                       2739, 10899, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23795, 3, 2739,
                                                                       2742, 10905, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23805, 3, 2742,
                                                                       2745, 10911, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23815, 3, 2745,
                                                                       2748, 10917, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23825, 3, 2748,
                                                                       2751, 10923, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23835, 3, 2751,
                                                                       2754, 10929, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23845, 3, 2754,
                                                                       2757, 10935, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23855, 3, 2757,
                                                                       2760, 10941, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23865, 3, 2760,
                                                                       2763, 10947, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23875, 3, 2763,
                                                                       2766, 10953, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23885, 3, 2766,
                                                                       2769, 10959, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 23895, 0, 3,
                                                                       23745, 10875, 23755, 2775,
                                                                       2784, 10965, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 23925, 0, 3,
                                                                       23755, 10881, 23765, 2784,
                                                                       2793, 10983, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 23955, 0, 3,
                                                                       23765, 10887, 23775, 2793,
                                                                       2802, 11001, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 23985, 0, 3,
                                                                       23775, 10893, 23785, 2802,
                                                                       2811, 11019, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24015, 0, 3,
                                                                       23785, 10899, 23795, 2811,
                                                                       2820, 11037, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24045, 0, 3,
                                                                       23795, 10905, 23805, 2820,
                                                                       2829, 11055, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24075, 0, 3,
                                                                       23805, 10911, 23815, 2829,
                                                                       2838, 11073, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24105, 0, 3,
                                                                       23815, 10917, 23825, 2838,
                                                                       2847, 11091, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24135, 0, 3,
                                                                       23825, 10923, 23835, 2847,
                                                                       2856, 11109, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24165, 0, 3,
                                                                       23835, 10929, 23845, 2856,
                                                                       2865, 11127, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24195, 0, 3,
                                                                       23845, 10935, 23855, 2865,
                                                                       2874, 11145, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24225, 0, 3,
                                                                       23855, 10941, 23865, 2874,
                                                                       2883, 11163, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24255, 0, 3,
                                                                       23865, 10947, 23875, 2883,
                                                                       2892, 11181, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24285, 0, 3,
                                                                       23875, 10953, 23885, 2892,
                                                                       2901, 11199, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24315, 0, 3,
                                                                       23895, 10965, 23925, 2919,
                                                                       2937, 11217, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24375, 0, 3,
                                                                       23925, 10983, 23955, 2937,
                                                                       2955, 11253, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24435, 0, 3,
                                                                       23955, 11001, 23985, 2955,
                                                                       2973, 11289, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24495, 0, 3,
                                                                       23985, 11019, 24015, 2973,
                                                                       2991, 11325, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24555, 0, 3,
                                                                       24015, 11037, 24045, 2991,
                                                                       3009, 11361, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24615, 0, 3,
                                                                       24045, 11055, 24075, 3009,
                                                                       3027, 11397, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24675, 0, 3,
                                                                       24075, 11073, 24105, 3027,
                                                                       3045, 11433, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24735, 0, 3,
                                                                       24105, 11091, 24135, 3045,
                                                                       3063, 11469, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24795, 0, 3,
                                                                       24135, 11109, 24165, 3063,
                                                                       3081, 11505, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24855, 0, 3,
                                                                       24165, 11127, 24195, 3081,
                                                                       3099, 11541, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24915, 0, 3,
                                                                       24195, 11145, 24225, 3099,
                                                                       3117, 11577, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24975, 0, 3,
                                                                       24225, 11163, 24255, 3117,
                                                                       3135, 11613, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 25035, 0, 3,
                                                                       24255, 11181, 24285, 3135,
                                                                       3153, 11649, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25095, 0, 3,
                                                                       24315, 11217, 24375, 3189,
                                                                       3219, 11685, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25195, 0, 3,
                                                                       24375, 11253, 24435, 3219,
                                                                       3249, 11745, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25295, 0, 3,
                                                                       24435, 11289, 24495, 3249,
                                                                       3279, 11805, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25395, 0, 3,
                                                                       24495, 11325, 24555, 3279,
                                                                       3309, 11865, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25495, 0, 3,
                                                                       24555, 11361, 24615, 3309,
                                                                       3339, 11925, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25595, 0, 3,
                                                                       24615, 11397, 24675, 3339,
                                                                       3369, 11985, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25695, 0, 3,
                                                                       24675, 11433, 24735, 3369,
                                                                       3399, 12045, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25795, 0, 3,
                                                                       24735, 11469, 24795, 3399,
                                                                       3429, 12105, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25895, 0, 3,
                                                                       24795, 11505, 24855, 3429,
                                                                       3459, 12165, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25995, 0, 3,
                                                                       24855, 11541, 24915, 3459,
                                                                       3489, 12225, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 26095, 0, 3,
                                                                       24915, 11577, 24975, 3489,
                                                                       3519, 12285, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 26195, 0, 3,
                                                                       24975, 11613, 25035, 3519,
                                                                       3549, 12345, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26295, 0, 3,
                                                                       25095, 11685, 25195, 3609,
                                                                       3654, 12405, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26445, 0, 3,
                                                                       25195, 11745, 25295, 3654,
                                                                       3699, 12495, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26595, 0, 3,
                                                                       25295, 11805, 25395, 3699,
                                                                       3744, 12585, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26745, 0, 3,
                                                                       25395, 11865, 25495, 3744,
                                                                       3789, 12675, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26895, 0, 3,
                                                                       25495, 11925, 25595, 3789,
                                                                       3834, 12765, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27045, 0, 3,
                                                                       25595, 11985, 25695, 3834,
                                                                       3879, 12855, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27195, 0, 3,
                                                                       25695, 12045, 25795, 3879,
                                                                       3924, 12945, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27345, 0, 3,
                                                                       25795, 12105, 25895, 3924,
                                                                       3969, 13035, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27495, 0, 3,
                                                                       25895, 12165, 25995, 3969,
                                                                       4014, 13125, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27645, 0, 3,
                                                                       25995, 12225, 26095, 4014,
                                                                       4059, 13215, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27795, 0, 3,
                                                                       26095, 12285, 26195, 4059,
                                                                       4104, 13305, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 27945, 0, 3,
                                                                       26295, 12405, 26445, 4194,
                                                                       4257, 13395, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 28155, 0, 3,
                                                                       26445, 12495, 26595, 4257,
                                                                       4320, 13521, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 28365, 0, 3,
                                                                       26595, 12585, 26745, 4320,
                                                                       4383, 13647, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 28575, 0, 3,
                                                                       26745, 12675, 26895, 4383,
                                                                       4446, 13773, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 28785, 0, 3,
                                                                       26895, 12765, 27045, 4446,
                                                                       4509, 13899, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 28995, 0, 3,
                                                                       27045, 12855, 27195, 4509,
                                                                       4572, 14025, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 29205, 0, 3,
                                                                       27195, 12945, 27345, 4572,
                                                                       4635, 14151, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 29415, 0, 3,
                                                                       27345, 13035, 27495, 4635,
                                                                       4698, 14277, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 29625, 0, 3,
                                                                       27495, 13125, 27645, 4698,
                                                                       4761, 14403, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 29835, 0, 3,
                                                                       27645, 13215, 27795, 4761,
                                                                       4824, 14529, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 30045, 0, 3,
                                                                       27945, 13395, 28155, 4950,
                                                                       5034, 14655, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 30325, 0, 3,
                                                                       28155, 13521, 28365, 5034,
                                                                       5118, 14823, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 30605, 0, 3,
                                                                       28365, 13647, 28575, 5118,
                                                                       5202, 14991, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 30885, 0, 3,
                                                                       28575, 13773, 28785, 5202,
                                                                       5286, 15159, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 31165, 0, 3,
                                                                       28785, 13899, 28995, 5286,
                                                                       5370, 15327, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 31445, 0, 3,
                                                                       28995, 14025, 29205, 5370,
                                                                       5454, 15495, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 31725, 0, 3,
                                                                       29205, 14151, 29415, 5454,
                                                                       5538, 15663, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 32005, 0, 3,
                                                                       29415, 14277, 29625, 5538,
                                                                       5622, 15831, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 32285, 0, 3,
                                                                       29625, 14403, 29835, 5622,
                                                                       5706, 15999, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 32565, 0, 3,
                                                                       30045, 14655, 30325, 5874,
                                                                       5982, 16167, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 32925, 0, 3,
                                                                       30325, 14823, 30605, 5982,
                                                                       6090, 16383, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 33285, 0, 3,
                                                                       30605, 14991, 30885, 6090,
                                                                       6198, 16599, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 33645, 0, 3,
                                                                       30885, 15159, 31165, 6198,
                                                                       6306, 16815, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 34005, 0, 3,
                                                                       31165, 15327, 31445, 6306,
                                                                       6414, 17031, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 34365, 0, 3,
                                                                       31445, 15495, 31725, 6414,
                                                                       6522, 17247, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 34725, 0, 3,
                                                                       31725, 15663, 32005, 6522,
                                                                       6630, 17463, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 35085, 0, 3,
                                                                       32005, 15831, 32285, 6630,
                                                                       6738, 17679, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 35445, 0, 3,
                                                                       32565, 16167, 32925, 6954,
                                                                       7089, 17895, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 35895, 0, 3,
                                                                       32925, 16383, 33285, 7089,
                                                                       7224, 18165, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 36345, 0, 3,
                                                                       33285, 16599, 33645, 7224,
                                                                       7359, 18435, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 36795, 0, 3,
                                                                       33645, 16815, 34005, 7359,
                                                                       7494, 18705, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 37245, 0, 3,
                                                                       34005, 17031, 34365, 7494,
                                                                       7629, 18975, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 37695, 0, 3,
                                                                       34365, 17247, 34725, 7629,
                                                                       7764, 19245, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 38145, 0, 3,
                                                                       34725, 17463, 35085, 7764,
                                                                       7899, 19515, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 38595, 0, 3,
                                                                       35445, 17895, 35895, 8169,
                                                                       8334, 19785, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 39145, 0, 3,
                                                                       35895, 18165, 36345, 8334,
                                                                       8499, 20115, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 39695, 0, 3,
                                                                       36345, 18435, 36795, 8499,
                                                                       8664, 20445, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 40245, 0, 3,
                                                                       36795, 18705, 37245, 8664,
                                                                       8829, 20775, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 40795, 0, 3,
                                                                       37245, 18975, 37695, 8829,
                                                                       8994, 21105, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 41345, 0, 3,
                                                                       37695, 19245, 38145, 8994,
                                                                       9159, 21435, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 41895, 0, 3,
                                                                       38595, 19785, 39145, 9489,
                                                                       9687, 21765, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 42555, 0, 3,
                                                                       39145, 20115, 39695, 9687,
                                                                       9885, 22161, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 43215, 0, 3,
                                                                       39695, 20445, 40245, 9885,
                                                                       10083, 22557, ncols,
                                                                       gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 43875, 0, 3,
                                                                       40245, 20775, 40795,
                                                                       10083, 10281, 22953,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 44535, 0, 3,
                                                                       40795, 21105, 41345,
                                                                       10281, 10479, 23349,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45195, 3, 10875,
                                                                       10881, 23765, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45210, 3, 10881,
                                                                       10887, 23775, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45225, 3, 10887,
                                                                       10893, 23785, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45240, 3, 10893,
                                                                       10899, 23795, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45255, 3, 10899,
                                                                       10905, 23805, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45270, 3, 10905,
                                                                       10911, 23815, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45285, 3, 10911,
                                                                       10917, 23825, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45300, 3, 10917,
                                                                       10923, 23835, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45315, 3, 10923,
                                                                       10929, 23845, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45330, 3, 10929,
                                                                       10935, 23855, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45345, 3, 10935,
                                                                       10941, 23865, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45360, 3, 10941,
                                                                       10947, 23875, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45375, 3, 10947,
                                                                       10953, 23885, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45390, 0, 3,
                                                                       45195, 23765, 45210,
                                                                       10965, 10983, 23955,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45435, 0, 3,
                                                                       45210, 23775, 45225,
                                                                       10983, 11001, 23985,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45480, 0, 3,
                                                                       45225, 23785, 45240,
                                                                       11001, 11019, 24015,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45525, 0, 3,
                                                                       45240, 23795, 45255,
                                                                       11019, 11037, 24045,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45570, 0, 3,
                                                                       45255, 23805, 45270,
                                                                       11037, 11055, 24075,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45615, 0, 3,
                                                                       45270, 23815, 45285,
                                                                       11055, 11073, 24105,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45660, 0, 3,
                                                                       45285, 23825, 45300,
                                                                       11073, 11091, 24135,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45705, 0, 3,
                                                                       45300, 23835, 45315,
                                                                       11091, 11109, 24165,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45750, 0, 3,
                                                                       45315, 23845, 45330,
                                                                       11109, 11127, 24195,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45795, 0, 3,
                                                                       45330, 23855, 45345,
                                                                       11127, 11145, 24225,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45840, 0, 3,
                                                                       45345, 23865, 45360,
                                                                       11145, 11163, 24255,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45885, 0, 3,
                                                                       45360, 23875, 45375,
                                                                       11163, 11181, 24285,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 45930, 0, 3,
                                                                       45390, 23955, 45435,
                                                                       11217, 11253, 24435,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46020, 0, 3,
                                                                       45435, 23985, 45480,
                                                                       11253, 11289, 24495,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46110, 0, 3,
                                                                       45480, 24015, 45525,
                                                                       11289, 11325, 24555,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46200, 0, 3,
                                                                       45525, 24045, 45570,
                                                                       11325, 11361, 24615,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46290, 0, 3,
                                                                       45570, 24075, 45615,
                                                                       11361, 11397, 24675,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46380, 0, 3,
                                                                       45615, 24105, 45660,
                                                                       11397, 11433, 24735,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46470, 0, 3,
                                                                       45660, 24135, 45705,
                                                                       11433, 11469, 24795,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46560, 0, 3,
                                                                       45705, 24165, 45750,
                                                                       11469, 11505, 24855,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46650, 0, 3,
                                                                       45750, 24195, 45795,
                                                                       11505, 11541, 24915,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46740, 0, 3,
                                                                       45795, 24225, 45840,
                                                                       11541, 11577, 24975,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46830, 0, 3,
                                                                       45840, 24255, 45885,
                                                                       11577, 11613, 25035,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 46920, 0, 3,
                                                                       45930, 24435, 46020,
                                                                       11685, 11745, 25295,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47070, 0, 3,
                                                                       46020, 24495, 46110,
                                                                       11745, 11805, 25395,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47220, 0, 3,
                                                                       46110, 24555, 46200,
                                                                       11805, 11865, 25495,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47370, 0, 3,
                                                                       46200, 24615, 46290,
                                                                       11865, 11925, 25595,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47520, 0, 3,
                                                                       46290, 24675, 46380,
                                                                       11925, 11985, 25695,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47670, 0, 3,
                                                                       46380, 24735, 46470,
                                                                       11985, 12045, 25795,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47820, 0, 3,
                                                                       46470, 24795, 46560,
                                                                       12045, 12105, 25895,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47970, 0, 3,
                                                                       46560, 24855, 46650,
                                                                       12105, 12165, 25995,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 48120, 0, 3,
                                                                       46650, 24915, 46740,
                                                                       12165, 12225, 26095,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 48270, 0, 3,
                                                                       46740, 24975, 46830,
                                                                       12225, 12285, 26195,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 48420, 0, 3,
                                                                       46920, 25295, 47070,
                                                                       12405, 12495, 26595,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 48645, 0, 3,
                                                                       47070, 25395, 47220,
                                                                       12495, 12585, 26745,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 48870, 0, 3,
                                                                       47220, 25495, 47370,
                                                                       12585, 12675, 26895,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 49095, 0, 3,
                                                                       47370, 25595, 47520,
                                                                       12675, 12765, 27045,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 49320, 0, 3,
                                                                       47520, 25695, 47670,
                                                                       12765, 12855, 27195,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 49545, 0, 3,
                                                                       47670, 25795, 47820,
                                                                       12855, 12945, 27345,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 49770, 0, 3,
                                                                       47820, 25895, 47970,
                                                                       12945, 13035, 27495,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 49995, 0, 3,
                                                                       47970, 25995, 48120,
                                                                       13035, 13125, 27645,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 50220, 0, 3,
                                                                       48120, 26095, 48270,
                                                                       13125, 13215, 27795,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 50445, 0, 3,
                                                                       48420, 26595, 48645,
                                                                       13395, 13521, 28365,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 50760, 0, 3,
                                                                       48645, 26745, 48870,
                                                                       13521, 13647, 28575,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 51075, 0, 3,
                                                                       48870, 26895, 49095,
                                                                       13647, 13773, 28785,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 51390, 0, 3,
                                                                       49095, 27045, 49320,
                                                                       13773, 13899, 28995,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 51705, 0, 3,
                                                                       49320, 27195, 49545,
                                                                       13899, 14025, 29205,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 52020, 0, 3,
                                                                       49545, 27345, 49770,
                                                                       14025, 14151, 29415,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 52335, 0, 3,
                                                                       49770, 27495, 49995,
                                                                       14151, 14277, 29625,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 52650, 0, 3,
                                                                       49995, 27645, 50220,
                                                                       14277, 14403, 29835,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 52965, 0, 3,
                                                                       50445, 28365, 50760,
                                                                       14655, 14823, 30605,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 53385, 0, 3,
                                                                       50760, 28575, 51075,
                                                                       14823, 14991, 30885,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 53805, 0, 3,
                                                                       51075, 28785, 51390,
                                                                       14991, 15159, 31165,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 54225, 0, 3,
                                                                       51390, 28995, 51705,
                                                                       15159, 15327, 31445,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 54645, 0, 3,
                                                                       51705, 29205, 52020,
                                                                       15327, 15495, 31725,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 55065, 0, 3,
                                                                       52020, 29415, 52335,
                                                                       15495, 15663, 32005,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 55485, 0, 3,
                                                                       52335, 29625, 52650,
                                                                       15663, 15831, 32285,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 55905, 0, 3,
                                                                       52965, 30605, 53385,
                                                                       16167, 16383, 33285,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 56445, 0, 3,
                                                                       53385, 30885, 53805,
                                                                       16383, 16599, 33645,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 56985, 0, 3,
                                                                       53805, 31165, 54225,
                                                                       16599, 16815, 34005,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 57525, 0, 3,
                                                                       54225, 31445, 54645,
                                                                       16815, 17031, 34365,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 58065, 0, 3,
                                                                       54645, 31725, 55065,
                                                                       17031, 17247, 34725,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 58605, 0, 3,
                                                                       55065, 32005, 55485,
                                                                       17247, 17463, 35085,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 59145, 0, 3,
                                                                       55905, 33285, 56445,
                                                                       17895, 18165, 36345,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 59820, 0, 3,
                                                                       56445, 33645, 56985,
                                                                       18165, 18435, 36795,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 60495, 0, 3,
                                                                       56985, 34005, 57525,
                                                                       18435, 18705, 37245,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 61170, 0, 3,
                                                                       57525, 34365, 58065,
                                                                       18705, 18975, 37695,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 61845, 0, 3,
                                                                       58065, 34725, 58605,
                                                                       18975, 19245, 38145,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 62520, 0, 3,
                                                                       59145, 36345, 59820,
                                                                       19785, 20115, 39695,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 63345, 0, 3,
                                                                       59820, 36795, 60495,
                                                                       20115, 20445, 40245,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 64170, 0, 3,
                                                                       60495, 37245, 61170,
                                                                       20445, 20775, 40795,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 64995, 0, 3,
                                                                       61170, 37695, 61845,
                                                                       20775, 21105, 41345,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 65820, 0, 3,
                                                                       62520, 39695, 63345,
                                                                       21765, 22161, 43215,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 66810, 0, 3,
                                                                       63345, 40245, 64170,
                                                                       22161, 22557, 43875,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 67800, 0, 3,
                                                                       64170, 40795, 64995,
                                                                       22557, 22953, 44535,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68790, 3, 23745,
                                                                       23755, 45195, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68811, 3, 23755,
                                                                       23765, 45210, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68832, 3, 23765,
                                                                       23775, 45225, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68853, 3, 23775,
                                                                       23785, 45240, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68874, 3, 23785,
                                                                       23795, 45255, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68895, 3, 23795,
                                                                       23805, 45270, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68916, 3, 23805,
                                                                       23815, 45285, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68937, 3, 23815,
                                                                       23825, 45300, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68958, 3, 23825,
                                                                       23835, 45315, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68979, 3, 23835,
                                                                       23845, 45330, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 69000, 3, 23845,
                                                                       23855, 45345, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 69021, 3, 23855,
                                                                       23865, 45360, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 69042, 3, 23865,
                                                                       23875, 45375, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69063, 0, 3,
                                                                       68790, 45195, 68811,
                                                                       23895, 23925, 45390,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69126, 0, 3,
                                                                       68811, 45210, 68832,
                                                                       23925, 23955, 45435,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69189, 0, 3,
                                                                       68832, 45225, 68853,
                                                                       23955, 23985, 45480,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69252, 0, 3,
                                                                       68853, 45240, 68874,
                                                                       23985, 24015, 45525,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69315, 0, 3,
                                                                       68874, 45255, 68895,
                                                                       24015, 24045, 45570,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69378, 0, 3,
                                                                       68895, 45270, 68916,
                                                                       24045, 24075, 45615,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69441, 0, 3,
                                                                       68916, 45285, 68937,
                                                                       24075, 24105, 45660,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69504, 0, 3,
                                                                       68937, 45300, 68958,
                                                                       24105, 24135, 45705,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69567, 0, 3,
                                                                       68958, 45315, 68979,
                                                                       24135, 24165, 45750,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69630, 0, 3,
                                                                       68979, 45330, 69000,
                                                                       24165, 24195, 45795,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69693, 0, 3,
                                                                       69000, 45345, 69021,
                                                                       24195, 24225, 45840,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 69756, 0, 3,
                                                                       69021, 45360, 69042,
                                                                       24225, 24255, 45885,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 69819, 0, 3,
                                                                       69063, 45390, 69126,
                                                                       24315, 24375, 45930,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 69945, 0, 3,
                                                                       69126, 45435, 69189,
                                                                       24375, 24435, 46020,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 70071, 0, 3,
                                                                       69189, 45480, 69252,
                                                                       24435, 24495, 46110,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 70197, 0, 3,
                                                                       69252, 45525, 69315,
                                                                       24495, 24555, 46200,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 70323, 0, 3,
                                                                       69315, 45570, 69378,
                                                                       24555, 24615, 46290,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 70449, 0, 3,
                                                                       69378, 45615, 69441,
                                                                       24615, 24675, 46380,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 70575, 0, 3,
                                                                       69441, 45660, 69504,
                                                                       24675, 24735, 46470,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 70701, 0, 3,
                                                                       69504, 45705, 69567,
                                                                       24735, 24795, 46560,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 70827, 0, 3,
                                                                       69567, 45750, 69630,
                                                                       24795, 24855, 46650,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 70953, 0, 3,
                                                                       69630, 45795, 69693,
                                                                       24855, 24915, 46740,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 71079, 0, 3,
                                                                       69693, 45840, 69756,
                                                                       24915, 24975, 46830,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 71205, 0, 3,
                                                                       69819, 45930, 69945,
                                                                       25095, 25195, 46920,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 71415, 0, 3,
                                                                       69945, 46020, 70071,
                                                                       25195, 25295, 47070,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 71625, 0, 3,
                                                                       70071, 46110, 70197,
                                                                       25295, 25395, 47220,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 71835, 0, 3,
                                                                       70197, 46200, 70323,
                                                                       25395, 25495, 47370,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 72045, 0, 3,
                                                                       70323, 46290, 70449,
                                                                       25495, 25595, 47520,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 72255, 0, 3,
                                                                       70449, 46380, 70575,
                                                                       25595, 25695, 47670,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 72465, 0, 3,
                                                                       70575, 46470, 70701,
                                                                       25695, 25795, 47820,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 72675, 0, 3,
                                                                       70701, 46560, 70827,
                                                                       25795, 25895, 47970,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 72885, 0, 3,
                                                                       70827, 46650, 70953,
                                                                       25895, 25995, 48120,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 73095, 0, 3,
                                                                       70953, 46740, 71079,
                                                                       25995, 26095, 48270,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 73305, 0, 3,
                                                                       71205, 46920, 71415,
                                                                       26295, 26445, 48420,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 73620, 0, 3,
                                                                       71415, 47070, 71625,
                                                                       26445, 26595, 48645,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 73935, 0, 3,
                                                                       71625, 47220, 71835,
                                                                       26595, 26745, 48870,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 74250, 0, 3,
                                                                       71835, 47370, 72045,
                                                                       26745, 26895, 49095,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 74565, 0, 3,
                                                                       72045, 47520, 72255,
                                                                       26895, 27045, 49320,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 74880, 0, 3,
                                                                       72255, 47670, 72465,
                                                                       27045, 27195, 49545,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 75195, 0, 3,
                                                                       72465, 47820, 72675,
                                                                       27195, 27345, 49770,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 75510, 0, 3,
                                                                       72675, 47970, 72885,
                                                                       27345, 27495, 49995,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 75825, 0, 3,
                                                                       72885, 48120, 73095,
                                                                       27495, 27645, 50220,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 76140, 0, 3,
                                                                       73305, 48420, 73620,
                                                                       27945, 28155, 50445,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 76581, 0, 3,
                                                                       73620, 48645, 73935,
                                                                       28155, 28365, 50760,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 77022, 0, 3,
                                                                       73935, 48870, 74250,
                                                                       28365, 28575, 51075,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 77463, 0, 3,
                                                                       74250, 49095, 74565,
                                                                       28575, 28785, 51390,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 77904, 0, 3,
                                                                       74565, 49320, 74880,
                                                                       28785, 28995, 51705,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 78345, 0, 3,
                                                                       74880, 49545, 75195,
                                                                       28995, 29205, 52020,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 78786, 0, 3,
                                                                       75195, 49770, 75510,
                                                                       29205, 29415, 52335,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 79227, 0, 3,
                                                                       75510, 49995, 75825,
                                                                       29415, 29625, 52650,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 79668, 0, 3,
                                                                       76140, 50445, 76581,
                                                                       30045, 30325, 52965,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 80256, 0, 3,
                                                                       76581, 50760, 77022,
                                                                       30325, 30605, 53385,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 80844, 0, 3,
                                                                       77022, 51075, 77463,
                                                                       30605, 30885, 53805,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 81432, 0, 3,
                                                                       77463, 51390, 77904,
                                                                       30885, 31165, 54225,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 82020, 0, 3,
                                                                       77904, 51705, 78345,
                                                                       31165, 31445, 54645,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 82608, 0, 3,
                                                                       78345, 52020, 78786,
                                                                       31445, 31725, 55065,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 83196, 0, 3,
                                                                       78786, 52335, 79227,
                                                                       31725, 32005, 55485,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 83784, 0, 3,
                                                                       79668, 52965, 80256,
                                                                       32565, 32925, 55905,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 84540, 0, 3,
                                                                       80256, 53385, 80844,
                                                                       32925, 33285, 56445,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 85296, 0, 3,
                                                                       80844, 53805, 81432,
                                                                       33285, 33645, 56985,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 86052, 0, 3,
                                                                       81432, 54225, 82020,
                                                                       33645, 34005, 57525,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 86808, 0, 3,
                                                                       82020, 54645, 82608,
                                                                       34005, 34365, 58065,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 87564, 0, 3,
                                                                       82608, 55065, 83196,
                                                                       34365, 34725, 58605,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 88320, 0, 3,
                                                                       83784, 55905, 84540,
                                                                       35445, 35895, 59145,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 89265, 0, 3,
                                                                       84540, 56445, 85296,
                                                                       35895, 36345, 59820,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 90210, 0, 3,
                                                                       85296, 56985, 86052,
                                                                       36345, 36795, 60495,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 91155, 0, 3,
                                                                       86052, 57525, 86808,
                                                                       36795, 37245, 61170,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 92100, 0, 3,
                                                                       86808, 58065, 87564,
                                                                       37245, 37695, 61845,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 93045, 0, 3,
                                                                       88320, 59145, 89265,
                                                                       38595, 39145, 62520,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 94200, 0, 3,
                                                                       89265, 59820, 90210,
                                                                       39145, 39695, 63345,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 95355, 0, 3,
                                                                       90210, 60495, 91155,
                                                                       39695, 40245, 64170,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 96510, 0, 3,
                                                                       91155, 61170, 92100,
                                                                       40245, 40795, 64995,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 97665, 0, 3,
                                                                       93045, 62520, 94200,
                                                                       41895, 42555, 65820,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 99051, 0, 3,
                                                                       94200, 63345, 95355,
                                                                       42555, 43215, 66810,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 100437, 0, 3,
                                                                       95355, 64170, 96510,
                                                                       43215, 43875, 67800,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101823, 3, 45195,
                                                                       45210, 68832, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101851, 3, 45210,
                                                                       45225, 68853, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101879, 3, 45225,
                                                                       45240, 68874, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101907, 3, 45240,
                                                                       45255, 68895, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101935, 3, 45255,
                                                                       45270, 68916, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101963, 3, 45270,
                                                                       45285, 68937, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101991, 3, 45285,
                                                                       45300, 68958, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 102019, 3, 45300,
                                                                       45315, 68979, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 102047, 3, 45315,
                                                                       45330, 69000, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 102075, 3, 45330,
                                                                       45345, 69021, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 102103, 3, 45345,
                                                                       45360, 69042, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 102131, 0, 3,
                                                                       101823, 68832, 101851,
                                                                       45390, 45435, 69189,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 102215, 0, 3,
                                                                       101851, 68853, 101879,
                                                                       45435, 45480, 69252,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 102299, 0, 3,
                                                                       101879, 68874, 101907,
                                                                       45480, 45525, 69315,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 102383, 0, 3,
                                                                       101907, 68895, 101935,
                                                                       45525, 45570, 69378,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 102467, 0, 3,
                                                                       101935, 68916, 101963,
                                                                       45570, 45615, 69441,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 102551, 0, 3,
                                                                       101963, 68937, 101991,
                                                                       45615, 45660, 69504,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 102635, 0, 3,
                                                                       101991, 68958, 102019,
                                                                       45660, 45705, 69567,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 102719, 0, 3,
                                                                       102019, 68979, 102047,
                                                                       45705, 45750, 69630,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 102803, 0, 3,
                                                                       102047, 69000, 102075,
                                                                       45750, 45795, 69693,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 102887, 0, 3,
                                                                       102075, 69021, 102103,
                                                                       45795, 45840, 69756,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 102971, 0, 3,
                                                                       102131, 69189, 102215,
                                                                       45930, 46020, 70071,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 103139, 0, 3,
                                                                       102215, 69252, 102299,
                                                                       46020, 46110, 70197,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 103307, 0, 3,
                                                                       102299, 69315, 102383,
                                                                       46110, 46200, 70323,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 103475, 0, 3,
                                                                       102383, 69378, 102467,
                                                                       46200, 46290, 70449,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 103643, 0, 3,
                                                                       102467, 69441, 102551,
                                                                       46290, 46380, 70575,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 103811, 0, 3,
                                                                       102551, 69504, 102635,
                                                                       46380, 46470, 70701,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 103979, 0, 3,
                                                                       102635, 69567, 102719,
                                                                       46470, 46560, 70827,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 104147, 0, 3,
                                                                       102719, 69630, 102803,
                                                                       46560, 46650, 70953,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 104315, 0, 3,
                                                                       102803, 69693, 102887,
                                                                       46650, 46740, 71079,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 104483, 0, 3,
                                                                       102971, 70071, 103139,
                                                                       46920, 47070, 71625,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 104763, 0, 3,
                                                                       103139, 70197, 103307,
                                                                       47070, 47220, 71835,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 105043, 0, 3,
                                                                       103307, 70323, 103475,
                                                                       47220, 47370, 72045,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 105323, 0, 3,
                                                                       103475, 70449, 103643,
                                                                       47370, 47520, 72255,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 105603, 0, 3,
                                                                       103643, 70575, 103811,
                                                                       47520, 47670, 72465,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 105883, 0, 3,
                                                                       103811, 70701, 103979,
                                                                       47670, 47820, 72675,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 106163, 0, 3,
                                                                       103979, 70827, 104147,
                                                                       47820, 47970, 72885,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 106443, 0, 3,
                                                                       104147, 70953, 104315,
                                                                       47970, 48120, 73095,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 106723, 0, 3,
                                                                       104483, 71625, 104763,
                                                                       48420, 48645, 73935,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 107143, 0, 3,
                                                                       104763, 71835, 105043,
                                                                       48645, 48870, 74250,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 107563, 0, 3,
                                                                       105043, 72045, 105323,
                                                                       48870, 49095, 74565,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 107983, 0, 3,
                                                                       105323, 72255, 105603,
                                                                       49095, 49320, 74880,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 108403, 0, 3,
                                                                       105603, 72465, 105883,
                                                                       49320, 49545, 75195,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 108823, 0, 3,
                                                                       105883, 72675, 106163,
                                                                       49545, 49770, 75510,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 109243, 0, 3,
                                                                       106163, 72885, 106443,
                                                                       49770, 49995, 75825,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 109663, 0, 3,
                                                                       106723, 73935, 107143,
                                                                       50445, 50760, 77022,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 110251, 0, 3,
                                                                       107143, 74250, 107563,
                                                                       50760, 51075, 77463,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 110839, 0, 3,
                                                                       107563, 74565, 107983,
                                                                       51075, 51390, 77904,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 111427, 0, 3,
                                                                       107983, 74880, 108403,
                                                                       51390, 51705, 78345,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 112015, 0, 3,
                                                                       108403, 75195, 108823,
                                                                       51705, 52020, 78786,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 112603, 0, 3,
                                                                       108823, 75510, 109243,
                                                                       52020, 52335, 79227,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 113191, 0, 3,
                                                                       109663, 77022, 110251,
                                                                       52965, 53385, 80844,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 113975, 0, 3,
                                                                       110251, 77463, 110839,
                                                                       53385, 53805, 81432,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 114759, 0, 3,
                                                                       110839, 77904, 111427,
                                                                       53805, 54225, 82020,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 115543, 0, 3,
                                                                       111427, 78345, 112015,
                                                                       54225, 54645, 82608,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 116327, 0, 3,
                                                                       112015, 78786, 112603,
                                                                       54645, 55065, 83196,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 117111, 0, 3,
                                                                       113191, 80844, 113975,
                                                                       55905, 56445, 85296,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 118119, 0, 3,
                                                                       113975, 81432, 114759,
                                                                       56445, 56985, 86052,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 119127, 0, 3,
                                                                       114759, 82020, 115543,
                                                                       56985, 57525, 86808,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 120135, 0, 3,
                                                                       115543, 82608, 116327,
                                                                       57525, 58065, 87564,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 121143, 0, 3,
                                                                       117111, 85296, 118119,
                                                                       59145, 59820, 90210,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 122403, 0, 3,
                                                                       118119, 86052, 119127,
                                                                       59820, 60495, 91155,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 123663, 0, 3,
                                                                       119127, 86808, 120135,
                                                                       60495, 61170, 92100,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 124923, 0, 3,
                                                                       121143, 90210, 122403,
                                                                       62520, 63345, 95355,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 126463, 0, 3,
                                                                       122403, 91155, 123663,
                                                                       63345, 64170, 96510,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 128003, 0, 3,
                                                                       124923, 95355, 126463,
                                                                       65820, 66810, 100437,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 129851, 3, 68790,
                                                                       68811, 101823, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 129887, 3, 68811,
                                                                       68832, 101851, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 129923, 3, 68832,
                                                                       68853, 101879, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 129959, 3, 68853,
                                                                       68874, 101907, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 129995, 3, 68874,
                                                                       68895, 101935, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130031, 3, 68895,
                                                                       68916, 101963, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130067, 3, 68916,
                                                                       68937, 101991, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130103, 3, 68937,
                                                                       68958, 102019, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130139, 3, 68958,
                                                                       68979, 102047, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130175, 3, 68979,
                                                                       69000, 102075, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130211, 3, 69000,
                                                                       69021, 102103, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 130247, 0, 3,
                                                                       129851, 101823, 129887,
                                                                       69063, 69126, 102131,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 130355, 0, 3,
                                                                       129887, 101851, 129923,
                                                                       69126, 69189, 102215,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 130463, 0, 3,
                                                                       129923, 101879, 129959,
                                                                       69189, 69252, 102299,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 130571, 0, 3,
                                                                       129959, 101907, 129995,
                                                                       69252, 69315, 102383,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 130679, 0, 3,
                                                                       129995, 101935, 130031,
                                                                       69315, 69378, 102467,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 130787, 0, 3,
                                                                       130031, 101963, 130067,
                                                                       69378, 69441, 102551,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 130895, 0, 3,
                                                                       130067, 101991, 130103,
                                                                       69441, 69504, 102635,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 131003, 0, 3,
                                                                       130103, 102019, 130139,
                                                                       69504, 69567, 102719,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 131111, 0, 3,
                                                                       130139, 102047, 130175,
                                                                       69567, 69630, 102803,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 131219, 0, 3,
                                                                       130175, 102075, 130211,
                                                                       69630, 69693, 102887,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 131327, 0, 3,
                                                                       130247, 102131, 130355,
                                                                       69819, 69945, 102971,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 131543, 0, 3,
                                                                       130355, 102215, 130463,
                                                                       69945, 70071, 103139,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 131759, 0, 3,
                                                                       130463, 102299, 130571,
                                                                       70071, 70197, 103307,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 131975, 0, 3,
                                                                       130571, 102383, 130679,
                                                                       70197, 70323, 103475,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 132191, 0, 3,
                                                                       130679, 102467, 130787,
                                                                       70323, 70449, 103643,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 132407, 0, 3,
                                                                       130787, 102551, 130895,
                                                                       70449, 70575, 103811,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 132623, 0, 3,
                                                                       130895, 102635, 131003,
                                                                       70575, 70701, 103979,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 132839, 0, 3,
                                                                       131003, 102719, 131111,
                                                                       70701, 70827, 104147,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 133055, 0, 3,
                                                                       131111, 102803, 131219,
                                                                       70827, 70953, 104315,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 133271, 0, 3,
                                                                       131327, 102971, 131543,
                                                                       71205, 71415, 104483,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 133631, 0, 3,
                                                                       131543, 103139, 131759,
                                                                       71415, 71625, 104763,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 133991, 0, 3,
                                                                       131759, 103307, 131975,
                                                                       71625, 71835, 105043,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 134351, 0, 3,
                                                                       131975, 103475, 132191,
                                                                       71835, 72045, 105323,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 134711, 0, 3,
                                                                       132191, 103643, 132407,
                                                                       72045, 72255, 105603,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 135071, 0, 3,
                                                                       132407, 103811, 132623,
                                                                       72255, 72465, 105883,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 135431, 0, 3,
                                                                       132623, 103979, 132839,
                                                                       72465, 72675, 106163,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 135791, 0, 3,
                                                                       132839, 104147, 133055,
                                                                       72675, 72885, 106443,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 136151, 0, 3,
                                                                       133271, 104483, 133631,
                                                                       73305, 73620, 106723,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 136691, 0, 3,
                                                                       133631, 104763, 133991,
                                                                       73620, 73935, 107143,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 137231, 0, 3,
                                                                       133991, 105043, 134351,
                                                                       73935, 74250, 107563,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 137771, 0, 3,
                                                                       134351, 105323, 134711,
                                                                       74250, 74565, 107983,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 138311, 0, 3,
                                                                       134711, 105603, 135071,
                                                                       74565, 74880, 108403,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 138851, 0, 3,
                                                                       135071, 105883, 135431,
                                                                       74880, 75195, 108823,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 139391, 0, 3,
                                                                       135431, 106163, 135791,
                                                                       75195, 75510, 109243,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 139931, 0, 3,
                                                                       136151, 106723, 136691,
                                                                       76140, 76581, 109663,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 140687, 0, 3,
                                                                       136691, 107143, 137231,
                                                                       76581, 77022, 110251,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 141443, 0, 3,
                                                                       137231, 107563, 137771,
                                                                       77022, 77463, 110839,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 142199, 0, 3,
                                                                       137771, 107983, 138311,
                                                                       77463, 77904, 111427,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 142955, 0, 3,
                                                                       138311, 108403, 138851,
                                                                       77904, 78345, 112015,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 143711, 0, 3,
                                                                       138851, 108823, 139391,
                                                                       78345, 78786, 112603,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 144467, 0, 3,
                                                                       139931, 109663, 140687,
                                                                       79668, 80256, 113191,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 145475, 0, 3,
                                                                       140687, 110251, 141443,
                                                                       80256, 80844, 113975,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 146483, 0, 3,
                                                                       141443, 110839, 142199,
                                                                       80844, 81432, 114759,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 147491, 0, 3,
                                                                       142199, 111427, 142955,
                                                                       81432, 82020, 115543,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 148499, 0, 3,
                                                                       142955, 112015, 143711,
                                                                       82020, 82608, 116327,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 149507, 0, 3,
                                                                       144467, 113191, 145475,
                                                                       83784, 84540, 117111,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 150803, 0, 3,
                                                                       145475, 113975, 146483,
                                                                       84540, 85296, 118119,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 152099, 0, 3,
                                                                       146483, 114759, 147491,
                                                                       85296, 86052, 119127,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 153395, 0, 3,
                                                                       147491, 115543, 148499,
                                                                       86052, 86808, 120135,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 154691, 0, 3,
                                                                       149507, 117111, 150803,
                                                                       88320, 89265, 121143,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 156311, 0, 3,
                                                                       150803, 118119, 152099,
                                                                       89265, 90210, 122403,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 157931, 0, 3,
                                                                       152099, 119127, 153395,
                                                                       90210, 91155, 123663,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 159551, 0, 3,
                                                                       154691, 121143, 156311,
                                                                       93045, 94200, 124923,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 161531, 0, 3,
                                                                       156311, 122403, 157931,
                                                                       94200, 95355, 126463,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 163511, 0, 3,
                                                                       159551, 124923, 161531,
                                                                       97665, 99051, 128003,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 165887, 144467, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 167315, 149507, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 169151, 154691, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 171446, 159551, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 174251, 163511, 2376, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 166895, 165887, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 168611, 167315, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 170771, 169151, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 173426, 171446, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 176627, 174251, 66, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 177617, 166895, 168611, 15,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 178877, 168611, 170771, 15,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 180497, 170771, 173426, 15,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 182522, 173426, 176627, 15,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 184997, 177617, 178877, 15,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 187517, 178877, 180497, 15,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 190757, 180497, 182522, 15,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 194807, 184997, 187517, 15,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 199007, 187517, 190757, 15,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 204407, 194807, 199007, 15,
                                             nmax);

        simdtrf::transform_g_inner(buffer, 210707, 204407, 28, 15, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 210707, 135, nmax);
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
