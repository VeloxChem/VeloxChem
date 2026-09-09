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


#include "SimdThreeCenterElectronRepulsionRecIIH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIF.hpp"
#include "SimdTransferIG.hpp"
#include "SimdTransferIH.hpp"
#include "SimdTransferII.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKF.hpp"
#include "SimdTransferKG.hpp"
#include "SimdTransferKH.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLD.hpp"
#include "SimdTransferLF.hpp"
#include "SimdTransferLG.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransferMD.hpp"
#include "SimdTransferMF.hpp"
#include "SimdTransferMP.hpp"
#include "SimdTransferND.hpp"
#include "SimdTransferNP.hpp"
#include "SimdTransferOP.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_iih_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_iih_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 225316, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1859 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 225316, 123767, 11767, dimensions);

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

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2724, 0, 3, 1822,
                                                                       1877, 2262, 2328, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2802, 0, 3, 1877,
                                                                       1932, 2328, 2394, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2880, 0, 3, 1932,
                                                                       1987, 2394, 2460, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2958, 0, 3, 1987,
                                                                       2042, 2460, 2526, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3036, 0, 3, 2042,
                                                                       2097, 2526, 2592, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3114, 0, 3, 2097,
                                                                       2152, 2592, 2658, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 3192, 0, 3, 2262,
                                                                       2328, 2724, 2802, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 3283, 0, 3, 2328,
                                                                       2394, 2802, 2880, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 3374, 0, 3, 2394,
                                                                       2460, 2880, 2958, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 3465, 0, 3, 2460,
                                                                       2526, 2958, 3036, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 3556, 0, 3, 2526,
                                                                       2592, 3036, 3114, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3647, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3650, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3653, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3656, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3659, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3662, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3665, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3668, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3671, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3674, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3677, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3680, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3683, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3686, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3689, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3692, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3695, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3698, 3, 7, 24,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3707, 3, 8, 27,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3716, 3, 9, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3725, 3, 10, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3734, 3, 11, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3743, 3, 12, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3752, 3, 13, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3761, 3, 14, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3770, 3, 15, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3779, 3, 16, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3788, 3, 17, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3797, 3, 18, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3806, 3, 19, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3815, 3, 20, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3824, 3, 21, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3833, 3, 22, 69,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3842, 3, 24, 72,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3860, 3, 27, 78,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3878, 3, 30, 84,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3896, 3, 33, 90,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3914, 3, 36, 96,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3932, 3, 39, 102,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3950, 3, 42, 108,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3968, 3, 45, 114,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3986, 3, 48, 120,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4004, 3, 51, 126,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4022, 3, 54, 132,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4040, 3, 57, 138,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4058, 3, 60, 144,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4076, 3, 63, 150,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4094, 3, 66, 156,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4112, 3, 72, 162,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4142, 3, 78, 172,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4172, 3, 84, 182,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4202, 3, 90, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4232, 3, 96, 202,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4262, 3, 102, 212,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4292, 3, 108, 222,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4322, 3, 114, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4352, 3, 120, 242,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4382, 3, 126, 252,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4412, 3, 132, 262,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4442, 3, 138, 272,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4472, 3, 144, 282,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4502, 3, 150, 292,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4532, 3, 162, 302,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4577, 3, 172, 317,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4622, 3, 182, 332,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4667, 3, 192, 347,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4712, 3, 202, 362,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4757, 3, 212, 377,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4802, 3, 222, 392,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4847, 3, 232, 407,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4892, 3, 242, 422,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4937, 3, 252, 437,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4982, 3, 262, 452,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5027, 3, 272, 467,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5072, 3, 282, 482,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5117, 3, 302, 497,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5180, 3, 317, 518,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5243, 3, 332, 539,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5306, 3, 347, 560,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5369, 3, 362, 581,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5432, 3, 377, 602,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5495, 3, 392, 623,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5558, 3, 407, 644,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5621, 3, 422, 665,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5684, 3, 437, 686,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5747, 3, 452, 707,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5810, 3, 467, 728,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5873, 3, 497, 749,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5957, 3, 518, 777,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6041, 3, 539, 805,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6125, 3, 560, 833,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6209, 3, 581, 861,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6293, 3, 602, 889,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6377, 3, 623, 917,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6461, 3, 644, 945,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6545, 3, 665, 973,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6629, 3, 686,
                                                                       1001, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6713, 3, 707,
                                                                       1029, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6797, 3, 749,
                                                                       1057, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6905, 3, 777,
                                                                       1093, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7013, 3, 805,
                                                                       1129, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7121, 3, 833,
                                                                       1165, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7229, 3, 861,
                                                                       1201, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7337, 3, 889,
                                                                       1237, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7445, 3, 917,
                                                                       1273, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7553, 3, 945,
                                                                       1309, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7661, 3, 973,
                                                                       1345, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7769, 3, 1001,
                                                                       1381, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7877, 3, 1057,
                                                                       1417, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8012, 3, 1093,
                                                                       1462, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8147, 3, 1129,
                                                                       1507, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8282, 3, 1165,
                                                                       1552, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8417, 3, 1201,
                                                                       1597, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8552, 3, 1237,
                                                                       1642, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8687, 3, 1273,
                                                                       1687, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8822, 3, 1309,
                                                                       1732, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8957, 3, 1345,
                                                                       1777, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9092, 3, 1417,
                                                                       1822, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9257, 3, 1462,
                                                                       1877, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9422, 3, 1507,
                                                                       1932, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9587, 3, 1552,
                                                                       1987, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9752, 3, 1597,
                                                                       2042, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9917, 3, 1642,
                                                                       2097, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10082, 3, 1687,
                                                                       2152, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10247, 3, 1732,
                                                                       2207, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10412, 3, 1822,
                                                                       2262, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10610, 3, 1877,
                                                                       2328, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10808, 3, 1932,
                                                                       2394, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11006, 3, 1987,
                                                                       2460, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11204, 3, 2042,
                                                                       2526, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11402, 3, 2097,
                                                                       2592, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11600, 3, 2152,
                                                                       2658, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 11798, 3, 2262,
                                                                       2724, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 12032, 3, 2328,
                                                                       2802, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 12266, 3, 2394,
                                                                       2880, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 12500, 3, 2460,
                                                                       2958, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 12734, 3, 2526,
                                                                       3036, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 12968, 3, 2592,
                                                                       3114, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 13202, 3, 2724,
                                                                       3192, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 13475, 3, 2802,
                                                                       3283, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 13748, 3, 2880,
                                                                       3374, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 14021, 3, 2958,
                                                                       3465, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 14294, 3, 3036,
                                                                       3556, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14567, 3, 7, 8,
                                                                       3653, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14573, 3, 8, 9,
                                                                       3656, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14579, 3, 9, 10,
                                                                       3659, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14585, 3, 10, 11,
                                                                       3662, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14591, 3, 11, 12,
                                                                       3665, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14597, 3, 12, 13,
                                                                       3668, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14603, 3, 13, 14,
                                                                       3671, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14609, 3, 14, 15,
                                                                       3674, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14615, 3, 15, 16,
                                                                       3677, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14621, 3, 16, 17,
                                                                       3680, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14627, 3, 17, 18,
                                                                       3683, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14633, 3, 18, 19,
                                                                       3686, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14639, 3, 19, 20,
                                                                       3689, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14645, 3, 20, 21,
                                                                       3692, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14651, 3, 21, 22,
                                                                       3695, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14657, 0, 3,
                                                                       14567, 3653, 14573, 3716,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14675, 0, 3,
                                                                       14573, 3656, 14579, 3725,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14693, 0, 3,
                                                                       14579, 3659, 14585, 3734,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14711, 0, 3,
                                                                       14585, 3662, 14591, 3743,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14729, 0, 3,
                                                                       14591, 3665, 14597, 3752,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14747, 0, 3,
                                                                       14597, 3668, 14603, 3761,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14765, 0, 3,
                                                                       14603, 3671, 14609, 3770,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14783, 0, 3,
                                                                       14609, 3674, 14615, 3779,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14801, 0, 3,
                                                                       14615, 3677, 14621, 3788,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14819, 0, 3,
                                                                       14621, 3680, 14627, 3797,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14837, 0, 3,
                                                                       14627, 3683, 14633, 3806,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14855, 0, 3,
                                                                       14633, 3686, 14639, 3815,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14873, 0, 3,
                                                                       14639, 3689, 14645, 3824,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14891, 0, 3,
                                                                       14645, 3692, 14651, 3833,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14909, 0, 3,
                                                                       14657, 3716, 14675, 72,
                                                                       78, 3878, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14945, 0, 3,
                                                                       14675, 3725, 14693, 78,
                                                                       84, 3896, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14981, 0, 3,
                                                                       14693, 3734, 14711, 84,
                                                                       90, 3914, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15017, 0, 3,
                                                                       14711, 3743, 14729, 90,
                                                                       96, 3932, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15053, 0, 3,
                                                                       14729, 3752, 14747, 96,
                                                                       102, 3950, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15089, 0, 3,
                                                                       14747, 3761, 14765, 102,
                                                                       108, 3968, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15125, 0, 3,
                                                                       14765, 3770, 14783, 108,
                                                                       114, 3986, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15161, 0, 3,
                                                                       14783, 3779, 14801, 114,
                                                                       120, 4004, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15197, 0, 3,
                                                                       14801, 3788, 14819, 120,
                                                                       126, 4022, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15233, 0, 3,
                                                                       14819, 3797, 14837, 126,
                                                                       132, 4040, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15269, 0, 3,
                                                                       14837, 3806, 14855, 132,
                                                                       138, 4058, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15305, 0, 3,
                                                                       14855, 3815, 14873, 138,
                                                                       144, 4076, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15341, 0, 3,
                                                                       14873, 3824, 14891, 144,
                                                                       150, 4094, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15377, 0, 3,
                                                                       14909, 3878, 14945, 162,
                                                                       172, 4172, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15437, 0, 3,
                                                                       14945, 3896, 14981, 172,
                                                                       182, 4202, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15497, 0, 3,
                                                                       14981, 3914, 15017, 182,
                                                                       192, 4232, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15557, 0, 3,
                                                                       15017, 3932, 15053, 192,
                                                                       202, 4262, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15617, 0, 3,
                                                                       15053, 3950, 15089, 202,
                                                                       212, 4292, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15677, 0, 3,
                                                                       15089, 3968, 15125, 212,
                                                                       222, 4322, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15737, 0, 3,
                                                                       15125, 3986, 15161, 222,
                                                                       232, 4352, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15797, 0, 3,
                                                                       15161, 4004, 15197, 232,
                                                                       242, 4382, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15857, 0, 3,
                                                                       15197, 4022, 15233, 242,
                                                                       252, 4412, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15917, 0, 3,
                                                                       15233, 4040, 15269, 252,
                                                                       262, 4442, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15977, 0, 3,
                                                                       15269, 4058, 15305, 262,
                                                                       272, 4472, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16037, 0, 3,
                                                                       15305, 4076, 15341, 272,
                                                                       282, 4502, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16097, 0, 3,
                                                                       15377, 4172, 15437, 302,
                                                                       317, 4622, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16187, 0, 3,
                                                                       15437, 4202, 15497, 317,
                                                                       332, 4667, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16277, 0, 3,
                                                                       15497, 4232, 15557, 332,
                                                                       347, 4712, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16367, 0, 3,
                                                                       15557, 4262, 15617, 347,
                                                                       362, 4757, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16457, 0, 3,
                                                                       15617, 4292, 15677, 362,
                                                                       377, 4802, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16547, 0, 3,
                                                                       15677, 4322, 15737, 377,
                                                                       392, 4847, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16637, 0, 3,
                                                                       15737, 4352, 15797, 392,
                                                                       407, 4892, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16727, 0, 3,
                                                                       15797, 4382, 15857, 407,
                                                                       422, 4937, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16817, 0, 3,
                                                                       15857, 4412, 15917, 422,
                                                                       437, 4982, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16907, 0, 3,
                                                                       15917, 4442, 15977, 437,
                                                                       452, 5027, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16997, 0, 3,
                                                                       15977, 4472, 16037, 452,
                                                                       467, 5072, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17087, 0, 3,
                                                                       16097, 4622, 16187, 497,
                                                                       518, 5243, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17213, 0, 3,
                                                                       16187, 4667, 16277, 518,
                                                                       539, 5306, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17339, 0, 3,
                                                                       16277, 4712, 16367, 539,
                                                                       560, 5369, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17465, 0, 3,
                                                                       16367, 4757, 16457, 560,
                                                                       581, 5432, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17591, 0, 3,
                                                                       16457, 4802, 16547, 581,
                                                                       602, 5495, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17717, 0, 3,
                                                                       16547, 4847, 16637, 602,
                                                                       623, 5558, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17843, 0, 3,
                                                                       16637, 4892, 16727, 623,
                                                                       644, 5621, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17969, 0, 3,
                                                                       16727, 4937, 16817, 644,
                                                                       665, 5684, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18095, 0, 3,
                                                                       16817, 4982, 16907, 665,
                                                                       686, 5747, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18221, 0, 3,
                                                                       16907, 5027, 16997, 686,
                                                                       707, 5810, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18347, 0, 3,
                                                                       17087, 5243, 17213, 749,
                                                                       777, 6041, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18515, 0, 3,
                                                                       17213, 5306, 17339, 777,
                                                                       805, 6125, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18683, 0, 3,
                                                                       17339, 5369, 17465, 805,
                                                                       833, 6209, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18851, 0, 3,
                                                                       17465, 5432, 17591, 833,
                                                                       861, 6293, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19019, 0, 3,
                                                                       17591, 5495, 17717, 861,
                                                                       889, 6377, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19187, 0, 3,
                                                                       17717, 5558, 17843, 889,
                                                                       917, 6461, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19355, 0, 3,
                                                                       17843, 5621, 17969, 917,
                                                                       945, 6545, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19523, 0, 3,
                                                                       17969, 5684, 18095, 945,
                                                                       973, 6629, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19691, 0, 3,
                                                                       18095, 5747, 18221, 973,
                                                                       1001, 6713, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 19859, 0, 3,
                                                                       18347, 6041, 18515, 1057,
                                                                       1093, 7013, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20075, 0, 3,
                                                                       18515, 6125, 18683, 1093,
                                                                       1129, 7121, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20291, 0, 3,
                                                                       18683, 6209, 18851, 1129,
                                                                       1165, 7229, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20507, 0, 3,
                                                                       18851, 6293, 19019, 1165,
                                                                       1201, 7337, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20723, 0, 3,
                                                                       19019, 6377, 19187, 1201,
                                                                       1237, 7445, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20939, 0, 3,
                                                                       19187, 6461, 19355, 1237,
                                                                       1273, 7553, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21155, 0, 3,
                                                                       19355, 6545, 19523, 1273,
                                                                       1309, 7661, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21371, 0, 3,
                                                                       19523, 6629, 19691, 1309,
                                                                       1345, 7769, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 21587, 0, 3,
                                                                       19859, 7013, 20075, 1417,
                                                                       1462, 8147, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 21857, 0, 3,
                                                                       20075, 7121, 20291, 1462,
                                                                       1507, 8282, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22127, 0, 3,
                                                                       20291, 7229, 20507, 1507,
                                                                       1552, 8417, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22397, 0, 3,
                                                                       20507, 7337, 20723, 1552,
                                                                       1597, 8552, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22667, 0, 3,
                                                                       20723, 7445, 20939, 1597,
                                                                       1642, 8687, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22937, 0, 3,
                                                                       20939, 7553, 21155, 1642,
                                                                       1687, 8822, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23207, 0, 3,
                                                                       21155, 7661, 21371, 1687,
                                                                       1732, 8957, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 23477, 0, 3,
                                                                       21587, 8147, 21857, 1822,
                                                                       1877, 9422, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 23807, 0, 3,
                                                                       21857, 8282, 22127, 1877,
                                                                       1932, 9587, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24137, 0, 3,
                                                                       22127, 8417, 22397, 1932,
                                                                       1987, 9752, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24467, 0, 3,
                                                                       22397, 8552, 22667, 1987,
                                                                       2042, 9917, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24797, 0, 3,
                                                                       22667, 8687, 22937, 2042,
                                                                       2097, 10082, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 25127, 0, 3,
                                                                       22937, 8822, 23207, 2097,
                                                                       2152, 10247, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 25457, 0, 3,
                                                                       23477, 9422, 23807, 2262,
                                                                       2328, 10808, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 25853, 0, 3,
                                                                       23807, 9587, 24137, 2328,
                                                                       2394, 11006, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 26249, 0, 3,
                                                                       24137, 9752, 24467, 2394,
                                                                       2460, 11204, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 26645, 0, 3,
                                                                       24467, 9917, 24797, 2460,
                                                                       2526, 11402, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 27041, 0, 3,
                                                                       24797, 10082, 25127, 2526,
                                                                       2592, 11600, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 27437, 0, 3,
                                                                       25457, 10808, 25853, 2724,
                                                                       2802, 12266, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 27905, 0, 3,
                                                                       25853, 11006, 26249, 2802,
                                                                       2880, 12500, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 28373, 0, 3,
                                                                       26249, 11204, 26645, 2880,
                                                                       2958, 12734, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 28841, 0, 3,
                                                                       26645, 11402, 27041, 2958,
                                                                       3036, 12968, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 29309, 0, 3,
                                                                       27437, 12266, 27905, 3192,
                                                                       3283, 13748, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 29855, 0, 3,
                                                                       27905, 12500, 28373, 3283,
                                                                       3374, 14021, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 30401, 0, 3,
                                                                       28373, 12734, 28841, 3374,
                                                                       3465, 14294, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30947, 3, 3647,
                                                                       3650, 14567, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30957, 3, 3650,
                                                                       3653, 14573, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30967, 3, 3653,
                                                                       3656, 14579, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30977, 3, 3656,
                                                                       3659, 14585, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30987, 3, 3659,
                                                                       3662, 14591, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 30997, 3, 3662,
                                                                       3665, 14597, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31007, 3, 3665,
                                                                       3668, 14603, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31017, 3, 3668,
                                                                       3671, 14609, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31027, 3, 3671,
                                                                       3674, 14615, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31037, 3, 3674,
                                                                       3677, 14621, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31047, 3, 3677,
                                                                       3680, 14627, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31057, 3, 3680,
                                                                       3683, 14633, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31067, 3, 3683,
                                                                       3686, 14639, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31077, 3, 3686,
                                                                       3689, 14645, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31087, 3, 3689,
                                                                       3692, 14651, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31097, 0, 3,
                                                                       30947, 14567, 30957, 3698,
                                                                       3707, 14657, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31127, 0, 3,
                                                                       30957, 14573, 30967, 3707,
                                                                       3716, 14675, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31157, 0, 3,
                                                                       30967, 14579, 30977, 3716,
                                                                       3725, 14693, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31187, 0, 3,
                                                                       30977, 14585, 30987, 3725,
                                                                       3734, 14711, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31217, 0, 3,
                                                                       30987, 14591, 30997, 3734,
                                                                       3743, 14729, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31247, 0, 3,
                                                                       30997, 14597, 31007, 3743,
                                                                       3752, 14747, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31277, 0, 3,
                                                                       31007, 14603, 31017, 3752,
                                                                       3761, 14765, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31307, 0, 3,
                                                                       31017, 14609, 31027, 3761,
                                                                       3770, 14783, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31337, 0, 3,
                                                                       31027, 14615, 31037, 3770,
                                                                       3779, 14801, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31367, 0, 3,
                                                                       31037, 14621, 31047, 3779,
                                                                       3788, 14819, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31397, 0, 3,
                                                                       31047, 14627, 31057, 3788,
                                                                       3797, 14837, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31427, 0, 3,
                                                                       31057, 14633, 31067, 3797,
                                                                       3806, 14855, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31457, 0, 3,
                                                                       31067, 14639, 31077, 3806,
                                                                       3815, 14873, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 31487, 0, 3,
                                                                       31077, 14645, 31087, 3815,
                                                                       3824, 14891, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31517, 0, 3,
                                                                       31097, 14657, 31127, 3842,
                                                                       3860, 14909, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31577, 0, 3,
                                                                       31127, 14675, 31157, 3860,
                                                                       3878, 14945, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31637, 0, 3,
                                                                       31157, 14693, 31187, 3878,
                                                                       3896, 14981, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31697, 0, 3,
                                                                       31187, 14711, 31217, 3896,
                                                                       3914, 15017, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31757, 0, 3,
                                                                       31217, 14729, 31247, 3914,
                                                                       3932, 15053, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31817, 0, 3,
                                                                       31247, 14747, 31277, 3932,
                                                                       3950, 15089, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31877, 0, 3,
                                                                       31277, 14765, 31307, 3950,
                                                                       3968, 15125, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31937, 0, 3,
                                                                       31307, 14783, 31337, 3968,
                                                                       3986, 15161, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 31997, 0, 3,
                                                                       31337, 14801, 31367, 3986,
                                                                       4004, 15197, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32057, 0, 3,
                                                                       31367, 14819, 31397, 4004,
                                                                       4022, 15233, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32117, 0, 3,
                                                                       31397, 14837, 31427, 4022,
                                                                       4040, 15269, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32177, 0, 3,
                                                                       31427, 14855, 31457, 4040,
                                                                       4058, 15305, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 32237, 0, 3,
                                                                       31457, 14873, 31487, 4058,
                                                                       4076, 15341, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32297, 0, 3,
                                                                       31517, 14909, 31577, 4112,
                                                                       4142, 15377, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32397, 0, 3,
                                                                       31577, 14945, 31637, 4142,
                                                                       4172, 15437, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32497, 0, 3,
                                                                       31637, 14981, 31697, 4172,
                                                                       4202, 15497, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32597, 0, 3,
                                                                       31697, 15017, 31757, 4202,
                                                                       4232, 15557, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32697, 0, 3,
                                                                       31757, 15053, 31817, 4232,
                                                                       4262, 15617, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32797, 0, 3,
                                                                       31817, 15089, 31877, 4262,
                                                                       4292, 15677, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32897, 0, 3,
                                                                       31877, 15125, 31937, 4292,
                                                                       4322, 15737, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 32997, 0, 3,
                                                                       31937, 15161, 31997, 4322,
                                                                       4352, 15797, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33097, 0, 3,
                                                                       31997, 15197, 32057, 4352,
                                                                       4382, 15857, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33197, 0, 3,
                                                                       32057, 15233, 32117, 4382,
                                                                       4412, 15917, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33297, 0, 3,
                                                                       32117, 15269, 32177, 4412,
                                                                       4442, 15977, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 33397, 0, 3,
                                                                       32177, 15305, 32237, 4442,
                                                                       4472, 16037, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 33497, 0, 3,
                                                                       32297, 15377, 32397, 4532,
                                                                       4577, 16097, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 33647, 0, 3,
                                                                       32397, 15437, 32497, 4577,
                                                                       4622, 16187, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 33797, 0, 3,
                                                                       32497, 15497, 32597, 4622,
                                                                       4667, 16277, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 33947, 0, 3,
                                                                       32597, 15557, 32697, 4667,
                                                                       4712, 16367, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34097, 0, 3,
                                                                       32697, 15617, 32797, 4712,
                                                                       4757, 16457, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34247, 0, 3,
                                                                       32797, 15677, 32897, 4757,
                                                                       4802, 16547, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34397, 0, 3,
                                                                       32897, 15737, 32997, 4802,
                                                                       4847, 16637, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34547, 0, 3,
                                                                       32997, 15797, 33097, 4847,
                                                                       4892, 16727, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34697, 0, 3,
                                                                       33097, 15857, 33197, 4892,
                                                                       4937, 16817, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34847, 0, 3,
                                                                       33197, 15917, 33297, 4937,
                                                                       4982, 16907, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 34997, 0, 3,
                                                                       33297, 15977, 33397, 4982,
                                                                       5027, 16997, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35147, 0, 3,
                                                                       33497, 16097, 33647, 5117,
                                                                       5180, 17087, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35357, 0, 3,
                                                                       33647, 16187, 33797, 5180,
                                                                       5243, 17213, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35567, 0, 3,
                                                                       33797, 16277, 33947, 5243,
                                                                       5306, 17339, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35777, 0, 3,
                                                                       33947, 16367, 34097, 5306,
                                                                       5369, 17465, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35987, 0, 3,
                                                                       34097, 16457, 34247, 5369,
                                                                       5432, 17591, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 36197, 0, 3,
                                                                       34247, 16547, 34397, 5432,
                                                                       5495, 17717, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 36407, 0, 3,
                                                                       34397, 16637, 34547, 5495,
                                                                       5558, 17843, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 36617, 0, 3,
                                                                       34547, 16727, 34697, 5558,
                                                                       5621, 17969, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 36827, 0, 3,
                                                                       34697, 16817, 34847, 5621,
                                                                       5684, 18095, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 37037, 0, 3,
                                                                       34847, 16907, 34997, 5684,
                                                                       5747, 18221, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 37247, 0, 3,
                                                                       35147, 17087, 35357, 5873,
                                                                       5957, 18347, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 37527, 0, 3,
                                                                       35357, 17213, 35567, 5957,
                                                                       6041, 18515, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 37807, 0, 3,
                                                                       35567, 17339, 35777, 6041,
                                                                       6125, 18683, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38087, 0, 3,
                                                                       35777, 17465, 35987, 6125,
                                                                       6209, 18851, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38367, 0, 3,
                                                                       35987, 17591, 36197, 6209,
                                                                       6293, 19019, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38647, 0, 3,
                                                                       36197, 17717, 36407, 6293,
                                                                       6377, 19187, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38927, 0, 3,
                                                                       36407, 17843, 36617, 6377,
                                                                       6461, 19355, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 39207, 0, 3,
                                                                       36617, 17969, 36827, 6461,
                                                                       6545, 19523, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 39487, 0, 3,
                                                                       36827, 18095, 37037, 6545,
                                                                       6629, 19691, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 39767, 0, 3,
                                                                       37247, 18347, 37527, 6797,
                                                                       6905, 19859, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 40127, 0, 3,
                                                                       37527, 18515, 37807, 6905,
                                                                       7013, 20075, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 40487, 0, 3,
                                                                       37807, 18683, 38087, 7013,
                                                                       7121, 20291, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 40847, 0, 3,
                                                                       38087, 18851, 38367, 7121,
                                                                       7229, 20507, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 41207, 0, 3,
                                                                       38367, 19019, 38647, 7229,
                                                                       7337, 20723, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 41567, 0, 3,
                                                                       38647, 19187, 38927, 7337,
                                                                       7445, 20939, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 41927, 0, 3,
                                                                       38927, 19355, 39207, 7445,
                                                                       7553, 21155, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 42287, 0, 3,
                                                                       39207, 19523, 39487, 7553,
                                                                       7661, 21371, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 42647, 0, 3,
                                                                       39767, 19859, 40127, 7877,
                                                                       8012, 21587, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 43097, 0, 3,
                                                                       40127, 20075, 40487, 8012,
                                                                       8147, 21857, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 43547, 0, 3,
                                                                       40487, 20291, 40847, 8147,
                                                                       8282, 22127, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 43997, 0, 3,
                                                                       40847, 20507, 41207, 8282,
                                                                       8417, 22397, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 44447, 0, 3,
                                                                       41207, 20723, 41567, 8417,
                                                                       8552, 22667, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 44897, 0, 3,
                                                                       41567, 20939, 41927, 8552,
                                                                       8687, 22937, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 45347, 0, 3,
                                                                       41927, 21155, 42287, 8687,
                                                                       8822, 23207, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 45797, 0, 3,
                                                                       42647, 21587, 43097, 9092,
                                                                       9257, 23477, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 46347, 0, 3,
                                                                       43097, 21857, 43547, 9257,
                                                                       9422, 23807, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 46897, 0, 3,
                                                                       43547, 22127, 43997, 9422,
                                                                       9587, 24137, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 47447, 0, 3,
                                                                       43997, 22397, 44447, 9587,
                                                                       9752, 24467, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 47997, 0, 3,
                                                                       44447, 22667, 44897, 9752,
                                                                       9917, 24797, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 48547, 0, 3,
                                                                       44897, 22937, 45347, 9917,
                                                                       10082, 25127, ncols,
                                                                       gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 49097, 0, 3,
                                                                       45797, 23477, 46347,
                                                                       10412, 10610, 25457,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 49757, 0, 3,
                                                                       46347, 23807, 46897,
                                                                       10610, 10808, 25853,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 50417, 0, 3,
                                                                       46897, 24137, 47447,
                                                                       10808, 11006, 26249,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 51077, 0, 3,
                                                                       47447, 24467, 47997,
                                                                       11006, 11204, 26645,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 51737, 0, 3,
                                                                       47997, 24797, 48547,
                                                                       11204, 11402, 27041,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 52397, 0, 3,
                                                                       49097, 25457, 49757,
                                                                       11798, 12032, 27437,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 53177, 0, 3,
                                                                       49757, 25853, 50417,
                                                                       12032, 12266, 27905,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 53957, 0, 3,
                                                                       50417, 26249, 51077,
                                                                       12266, 12500, 28373,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 54737, 0, 3,
                                                                       51077, 26645, 51737,
                                                                       12500, 12734, 28841,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 55517, 0, 3,
                                                                       52397, 27437, 53177,
                                                                       13202, 13475, 29309,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 56427, 0, 3,
                                                                       53177, 27905, 53957,
                                                                       13475, 13748, 29855,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 57337, 0, 3,
                                                                       53957, 28373, 54737,
                                                                       13748, 14021, 30401,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58247, 3, 14567,
                                                                       14573, 30967, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58262, 3, 14573,
                                                                       14579, 30977, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58277, 3, 14579,
                                                                       14585, 30987, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58292, 3, 14585,
                                                                       14591, 30997, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58307, 3, 14591,
                                                                       14597, 31007, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58322, 3, 14597,
                                                                       14603, 31017, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58337, 3, 14603,
                                                                       14609, 31027, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58352, 3, 14609,
                                                                       14615, 31037, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58367, 3, 14615,
                                                                       14621, 31047, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58382, 3, 14621,
                                                                       14627, 31057, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58397, 3, 14627,
                                                                       14633, 31067, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58412, 3, 14633,
                                                                       14639, 31077, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 58427, 3, 14639,
                                                                       14645, 31087, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58442, 0, 3,
                                                                       58247, 30967, 58262,
                                                                       14657, 14675, 31157,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58487, 0, 3,
                                                                       58262, 30977, 58277,
                                                                       14675, 14693, 31187,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58532, 0, 3,
                                                                       58277, 30987, 58292,
                                                                       14693, 14711, 31217,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58577, 0, 3,
                                                                       58292, 30997, 58307,
                                                                       14711, 14729, 31247,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58622, 0, 3,
                                                                       58307, 31007, 58322,
                                                                       14729, 14747, 31277,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58667, 0, 3,
                                                                       58322, 31017, 58337,
                                                                       14747, 14765, 31307,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58712, 0, 3,
                                                                       58337, 31027, 58352,
                                                                       14765, 14783, 31337,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58757, 0, 3,
                                                                       58352, 31037, 58367,
                                                                       14783, 14801, 31367,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58802, 0, 3,
                                                                       58367, 31047, 58382,
                                                                       14801, 14819, 31397,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58847, 0, 3,
                                                                       58382, 31057, 58397,
                                                                       14819, 14837, 31427,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58892, 0, 3,
                                                                       58397, 31067, 58412,
                                                                       14837, 14855, 31457,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 58937, 0, 3,
                                                                       58412, 31077, 58427,
                                                                       14855, 14873, 31487,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 58982, 0, 3,
                                                                       58442, 31157, 58487,
                                                                       14909, 14945, 31637,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 59072, 0, 3,
                                                                       58487, 31187, 58532,
                                                                       14945, 14981, 31697,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 59162, 0, 3,
                                                                       58532, 31217, 58577,
                                                                       14981, 15017, 31757,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 59252, 0, 3,
                                                                       58577, 31247, 58622,
                                                                       15017, 15053, 31817,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 59342, 0, 3,
                                                                       58622, 31277, 58667,
                                                                       15053, 15089, 31877,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 59432, 0, 3,
                                                                       58667, 31307, 58712,
                                                                       15089, 15125, 31937,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 59522, 0, 3,
                                                                       58712, 31337, 58757,
                                                                       15125, 15161, 31997,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 59612, 0, 3,
                                                                       58757, 31367, 58802,
                                                                       15161, 15197, 32057,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 59702, 0, 3,
                                                                       58802, 31397, 58847,
                                                                       15197, 15233, 32117,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 59792, 0, 3,
                                                                       58847, 31427, 58892,
                                                                       15233, 15269, 32177,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 59882, 0, 3,
                                                                       58892, 31457, 58937,
                                                                       15269, 15305, 32237,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 59972, 0, 3,
                                                                       58982, 31637, 59072,
                                                                       15377, 15437, 32497,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 60122, 0, 3,
                                                                       59072, 31697, 59162,
                                                                       15437, 15497, 32597,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 60272, 0, 3,
                                                                       59162, 31757, 59252,
                                                                       15497, 15557, 32697,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 60422, 0, 3,
                                                                       59252, 31817, 59342,
                                                                       15557, 15617, 32797,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 60572, 0, 3,
                                                                       59342, 31877, 59432,
                                                                       15617, 15677, 32897,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 60722, 0, 3,
                                                                       59432, 31937, 59522,
                                                                       15677, 15737, 32997,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 60872, 0, 3,
                                                                       59522, 31997, 59612,
                                                                       15737, 15797, 33097,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 61022, 0, 3,
                                                                       59612, 32057, 59702,
                                                                       15797, 15857, 33197,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 61172, 0, 3,
                                                                       59702, 32117, 59792,
                                                                       15857, 15917, 33297,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 61322, 0, 3,
                                                                       59792, 32177, 59882,
                                                                       15917, 15977, 33397,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 61472, 0, 3,
                                                                       59972, 32497, 60122,
                                                                       16097, 16187, 33797,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 61697, 0, 3,
                                                                       60122, 32597, 60272,
                                                                       16187, 16277, 33947,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 61922, 0, 3,
                                                                       60272, 32697, 60422,
                                                                       16277, 16367, 34097,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 62147, 0, 3,
                                                                       60422, 32797, 60572,
                                                                       16367, 16457, 34247,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 62372, 0, 3,
                                                                       60572, 32897, 60722,
                                                                       16457, 16547, 34397,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 62597, 0, 3,
                                                                       60722, 32997, 60872,
                                                                       16547, 16637, 34547,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 62822, 0, 3,
                                                                       60872, 33097, 61022,
                                                                       16637, 16727, 34697,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 63047, 0, 3,
                                                                       61022, 33197, 61172,
                                                                       16727, 16817, 34847,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 63272, 0, 3,
                                                                       61172, 33297, 61322,
                                                                       16817, 16907, 34997,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 63497, 0, 3,
                                                                       61472, 33797, 61697,
                                                                       17087, 17213, 35567,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 63812, 0, 3,
                                                                       61697, 33947, 61922,
                                                                       17213, 17339, 35777,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 64127, 0, 3,
                                                                       61922, 34097, 62147,
                                                                       17339, 17465, 35987,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 64442, 0, 3,
                                                                       62147, 34247, 62372,
                                                                       17465, 17591, 36197,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 64757, 0, 3,
                                                                       62372, 34397, 62597,
                                                                       17591, 17717, 36407,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 65072, 0, 3,
                                                                       62597, 34547, 62822,
                                                                       17717, 17843, 36617,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 65387, 0, 3,
                                                                       62822, 34697, 63047,
                                                                       17843, 17969, 36827,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 65702, 0, 3,
                                                                       63047, 34847, 63272,
                                                                       17969, 18095, 37037,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 66017, 0, 3,
                                                                       63497, 35567, 63812,
                                                                       18347, 18515, 37807,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 66437, 0, 3,
                                                                       63812, 35777, 64127,
                                                                       18515, 18683, 38087,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 66857, 0, 3,
                                                                       64127, 35987, 64442,
                                                                       18683, 18851, 38367,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 67277, 0, 3,
                                                                       64442, 36197, 64757,
                                                                       18851, 19019, 38647,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 67697, 0, 3,
                                                                       64757, 36407, 65072,
                                                                       19019, 19187, 38927,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 68117, 0, 3,
                                                                       65072, 36617, 65387,
                                                                       19187, 19355, 39207,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 68537, 0, 3,
                                                                       65387, 36827, 65702,
                                                                       19355, 19523, 39487,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 68957, 0, 3,
                                                                       66017, 37807, 66437,
                                                                       19859, 20075, 40487,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 69497, 0, 3,
                                                                       66437, 38087, 66857,
                                                                       20075, 20291, 40847,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 70037, 0, 3,
                                                                       66857, 38367, 67277,
                                                                       20291, 20507, 41207,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 70577, 0, 3,
                                                                       67277, 38647, 67697,
                                                                       20507, 20723, 41567,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 71117, 0, 3,
                                                                       67697, 38927, 68117,
                                                                       20723, 20939, 41927,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 71657, 0, 3,
                                                                       68117, 39207, 68537,
                                                                       20939, 21155, 42287,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 72197, 0, 3,
                                                                       68957, 40487, 69497,
                                                                       21587, 21857, 43547,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 72872, 0, 3,
                                                                       69497, 40847, 70037,
                                                                       21857, 22127, 43997,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 73547, 0, 3,
                                                                       70037, 41207, 70577,
                                                                       22127, 22397, 44447,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 74222, 0, 3,
                                                                       70577, 41567, 71117,
                                                                       22397, 22667, 44897,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 74897, 0, 3,
                                                                       71117, 41927, 71657,
                                                                       22667, 22937, 45347,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 75572, 0, 3,
                                                                       72197, 43547, 72872,
                                                                       23477, 23807, 46897,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 76397, 0, 3,
                                                                       72872, 43997, 73547,
                                                                       23807, 24137, 47447,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 77222, 0, 3,
                                                                       73547, 44447, 74222,
                                                                       24137, 24467, 47997,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 78047, 0, 3,
                                                                       74222, 44897, 74897,
                                                                       24467, 24797, 48547,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 78872, 0, 3,
                                                                       75572, 46897, 76397,
                                                                       25457, 25853, 50417,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 79862, 0, 3,
                                                                       76397, 47447, 77222,
                                                                       25853, 26249, 51077,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 80852, 0, 3,
                                                                       77222, 47997, 78047,
                                                                       26249, 26645, 51737,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 81842, 0, 3,
                                                                       78872, 50417, 79862,
                                                                       27437, 27905, 53957,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 83012, 0, 3,
                                                                       79862, 51077, 80852,
                                                                       27905, 28373, 54737,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 84182, 0, 3,
                                                                       81842, 53957, 83012,
                                                                       29309, 29855, 57337,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85547, 3, 30947,
                                                                       30957, 58247, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85568, 3, 30957,
                                                                       30967, 58262, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85589, 3, 30967,
                                                                       30977, 58277, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85610, 3, 30977,
                                                                       30987, 58292, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85631, 3, 30987,
                                                                       30997, 58307, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85652, 3, 30997,
                                                                       31007, 58322, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85673, 3, 31007,
                                                                       31017, 58337, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85694, 3, 31017,
                                                                       31027, 58352, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85715, 3, 31027,
                                                                       31037, 58367, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85736, 3, 31037,
                                                                       31047, 58382, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85757, 3, 31047,
                                                                       31057, 58397, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85778, 3, 31057,
                                                                       31067, 58412, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 85799, 3, 31067,
                                                                       31077, 58427, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 85820, 0, 3,
                                                                       85547, 58247, 85568,
                                                                       31097, 31127, 58442,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 85883, 0, 3,
                                                                       85568, 58262, 85589,
                                                                       31127, 31157, 58487,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 85946, 0, 3,
                                                                       85589, 58277, 85610,
                                                                       31157, 31187, 58532,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 86009, 0, 3,
                                                                       85610, 58292, 85631,
                                                                       31187, 31217, 58577,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 86072, 0, 3,
                                                                       85631, 58307, 85652,
                                                                       31217, 31247, 58622,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 86135, 0, 3,
                                                                       85652, 58322, 85673,
                                                                       31247, 31277, 58667,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 86198, 0, 3,
                                                                       85673, 58337, 85694,
                                                                       31277, 31307, 58712,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 86261, 0, 3,
                                                                       85694, 58352, 85715,
                                                                       31307, 31337, 58757,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 86324, 0, 3,
                                                                       85715, 58367, 85736,
                                                                       31337, 31367, 58802,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 86387, 0, 3,
                                                                       85736, 58382, 85757,
                                                                       31367, 31397, 58847,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 86450, 0, 3,
                                                                       85757, 58397, 85778,
                                                                       31397, 31427, 58892,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 86513, 0, 3,
                                                                       85778, 58412, 85799,
                                                                       31427, 31457, 58937,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 86576, 0, 3,
                                                                       85820, 58442, 85883,
                                                                       31517, 31577, 58982,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 86702, 0, 3,
                                                                       85883, 58487, 85946,
                                                                       31577, 31637, 59072,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 86828, 0, 3,
                                                                       85946, 58532, 86009,
                                                                       31637, 31697, 59162,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 86954, 0, 3,
                                                                       86009, 58577, 86072,
                                                                       31697, 31757, 59252,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 87080, 0, 3,
                                                                       86072, 58622, 86135,
                                                                       31757, 31817, 59342,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 87206, 0, 3,
                                                                       86135, 58667, 86198,
                                                                       31817, 31877, 59432,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 87332, 0, 3,
                                                                       86198, 58712, 86261,
                                                                       31877, 31937, 59522,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 87458, 0, 3,
                                                                       86261, 58757, 86324,
                                                                       31937, 31997, 59612,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 87584, 0, 3,
                                                                       86324, 58802, 86387,
                                                                       31997, 32057, 59702,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 87710, 0, 3,
                                                                       86387, 58847, 86450,
                                                                       32057, 32117, 59792,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 87836, 0, 3,
                                                                       86450, 58892, 86513,
                                                                       32117, 32177, 59882,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 87962, 0, 3,
                                                                       86576, 58982, 86702,
                                                                       32297, 32397, 59972,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 88172, 0, 3,
                                                                       86702, 59072, 86828,
                                                                       32397, 32497, 60122,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 88382, 0, 3,
                                                                       86828, 59162, 86954,
                                                                       32497, 32597, 60272,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 88592, 0, 3,
                                                                       86954, 59252, 87080,
                                                                       32597, 32697, 60422,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 88802, 0, 3,
                                                                       87080, 59342, 87206,
                                                                       32697, 32797, 60572,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 89012, 0, 3,
                                                                       87206, 59432, 87332,
                                                                       32797, 32897, 60722,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 89222, 0, 3,
                                                                       87332, 59522, 87458,
                                                                       32897, 32997, 60872,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 89432, 0, 3,
                                                                       87458, 59612, 87584,
                                                                       32997, 33097, 61022,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 89642, 0, 3,
                                                                       87584, 59702, 87710,
                                                                       33097, 33197, 61172,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 89852, 0, 3,
                                                                       87710, 59792, 87836,
                                                                       33197, 33297, 61322,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 90062, 0, 3,
                                                                       87962, 59972, 88172,
                                                                       33497, 33647, 61472,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 90377, 0, 3,
                                                                       88172, 60122, 88382,
                                                                       33647, 33797, 61697,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 90692, 0, 3,
                                                                       88382, 60272, 88592,
                                                                       33797, 33947, 61922,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 91007, 0, 3,
                                                                       88592, 60422, 88802,
                                                                       33947, 34097, 62147,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 91322, 0, 3,
                                                                       88802, 60572, 89012,
                                                                       34097, 34247, 62372,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 91637, 0, 3,
                                                                       89012, 60722, 89222,
                                                                       34247, 34397, 62597,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 91952, 0, 3,
                                                                       89222, 60872, 89432,
                                                                       34397, 34547, 62822,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 92267, 0, 3,
                                                                       89432, 61022, 89642,
                                                                       34547, 34697, 63047,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 92582, 0, 3,
                                                                       89642, 61172, 89852,
                                                                       34697, 34847, 63272,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 92897, 0, 3,
                                                                       90062, 61472, 90377,
                                                                       35147, 35357, 63497,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 93338, 0, 3,
                                                                       90377, 61697, 90692,
                                                                       35357, 35567, 63812,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 93779, 0, 3,
                                                                       90692, 61922, 91007,
                                                                       35567, 35777, 64127,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 94220, 0, 3,
                                                                       91007, 62147, 91322,
                                                                       35777, 35987, 64442,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 94661, 0, 3,
                                                                       91322, 62372, 91637,
                                                                       35987, 36197, 64757,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 95102, 0, 3,
                                                                       91637, 62597, 91952,
                                                                       36197, 36407, 65072,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 95543, 0, 3,
                                                                       91952, 62822, 92267,
                                                                       36407, 36617, 65387,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 95984, 0, 3,
                                                                       92267, 63047, 92582,
                                                                       36617, 36827, 65702,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 96425, 0, 3,
                                                                       92897, 63497, 93338,
                                                                       37247, 37527, 66017,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 97013, 0, 3,
                                                                       93338, 63812, 93779,
                                                                       37527, 37807, 66437,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 97601, 0, 3,
                                                                       93779, 64127, 94220,
                                                                       37807, 38087, 66857,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 98189, 0, 3,
                                                                       94220, 64442, 94661,
                                                                       38087, 38367, 67277,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 98777, 0, 3,
                                                                       94661, 64757, 95102,
                                                                       38367, 38647, 67697,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 99365, 0, 3,
                                                                       95102, 65072, 95543,
                                                                       38647, 38927, 68117,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 99953, 0, 3,
                                                                       95543, 65387, 95984,
                                                                       38927, 39207, 68537,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 100541, 0, 3,
                                                                       96425, 66017, 97013,
                                                                       39767, 40127, 68957,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 101297, 0, 3,
                                                                       97013, 66437, 97601,
                                                                       40127, 40487, 69497,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 102053, 0, 3,
                                                                       97601, 66857, 98189,
                                                                       40487, 40847, 70037,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 102809, 0, 3,
                                                                       98189, 67277, 98777,
                                                                       40847, 41207, 70577,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 103565, 0, 3,
                                                                       98777, 67697, 99365,
                                                                       41207, 41567, 71117,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 104321, 0, 3,
                                                                       99365, 68117, 99953,
                                                                       41567, 41927, 71657,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 105077, 0, 3,
                                                                       100541, 68957, 101297,
                                                                       42647, 43097, 72197,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 106022, 0, 3,
                                                                       101297, 69497, 102053,
                                                                       43097, 43547, 72872,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 106967, 0, 3,
                                                                       102053, 70037, 102809,
                                                                       43547, 43997, 73547,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 107912, 0, 3,
                                                                       102809, 70577, 103565,
                                                                       43997, 44447, 74222,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 108857, 0, 3,
                                                                       103565, 71117, 104321,
                                                                       44447, 44897, 74897,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 109802, 0, 3,
                                                                       105077, 72197, 106022,
                                                                       45797, 46347, 75572,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 110957, 0, 3,
                                                                       106022, 72872, 106967,
                                                                       46347, 46897, 76397,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 112112, 0, 3,
                                                                       106967, 73547, 107912,
                                                                       46897, 47447, 77222,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 113267, 0, 3,
                                                                       107912, 74222, 108857,
                                                                       47447, 47997, 78047,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 114422, 0, 3,
                                                                       109802, 75572, 110957,
                                                                       49097, 49757, 78872,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 115808, 0, 3,
                                                                       110957, 76397, 112112,
                                                                       49757, 50417, 79862,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 117194, 0, 3,
                                                                       112112, 77222, 113267,
                                                                       50417, 51077, 80852,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 118580, 0, 3,
                                                                       114422, 78872, 115808,
                                                                       52397, 53177, 81842,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 120218, 0, 3,
                                                                       115808, 79862, 117194,
                                                                       53177, 53957, 83012,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 121856, 0, 3,
                                                                       118580, 81842, 120218,
                                                                       55517, 56427, 84182,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 123767, 96425, 588, ncols);

                    simdfunc::contract_primitives(buffer, 124663, 100541, 756, ncols);

                    simdfunc::contract_primitives(buffer, 125815, 105077, 945, ncols);

                    simdfunc::contract_primitives(buffer, 127255, 109802, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 129015, 114422, 1386, ncols);

                    simdfunc::contract_primitives(buffer, 131127, 118580, 1638, ncols);

                    simdfunc::contract_primitives(buffer, 133623, 121856, 1911, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 124355, 123767, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 125419, 124663, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 126760, 125815, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 128410, 127255, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 130401, 129015, 66, 1, nmax);

        simdtrf::transform_h_inner(buffer, 132765, 131127, 78, 1, nmax);

        simdtrf::transform_h_inner(buffer, 135534, 133623, 91, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 136535, 124355, 125419, 11,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 137459, 125419, 126760, 11,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 138647, 126760, 128410, 11,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 140132, 128410, 130401, 11,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 141947, 130401, 132765, 11,
                                             nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 144125, 132765, 135534, 11,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 146699, 136535, 137459, 11,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 148547, 137459, 138647, 11,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 150923, 138647, 140132, 11,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 153893, 140132, 141947, 11,
                                             nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 157523, 141947, 144125, 11,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 161879, 146699, 148547, 11,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 164959, 148547, 150923, 11,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 168919, 150923, 153893, 11,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 173869, 153893, 157523, 11,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 179919, 161879, 164959, 11,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 184539, 164959, 168919, 11,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 190479, 168919, 173869, 11,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 197904, 179919, 184539, 11,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 204372, 184539, 190479, 11,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 212688, 197904, 204372, 11,
                                             nmax);

        simdtrf::transform_i_inner(buffer, 221312, 212688, 28, 11, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 221312, 143, nmax);
    }

    for (size_t m = 0; m < 1859; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
