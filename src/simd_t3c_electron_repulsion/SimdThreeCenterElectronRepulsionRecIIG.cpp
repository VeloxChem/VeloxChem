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


#include "SimdThreeCenterElectronRepulsionRecIIG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
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
#include "SimdTransformG.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_iig_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_iig_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 155932, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 155932, 73717, 8757, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 16,
                                                             ncols, fj, mu, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3647, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3650, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3653, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3656, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3659, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3662, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3665, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3668, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3671, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3674, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3677, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3680, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3683, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3686, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3689, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3692, 3, 9, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3701, 3, 10, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3710, 3, 11, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3719, 3, 12, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3728, 3, 13, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3737, 3, 14, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3746, 3, 15, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3755, 3, 16, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3764, 3, 17, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3773, 3, 18, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3782, 3, 19, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3791, 3, 20, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3800, 3, 21, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3809, 3, 22, 69,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3818, 3, 30, 84,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3836, 3, 33, 90,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3854, 3, 36, 96,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3872, 3, 39, 102,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3890, 3, 42, 108,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3908, 3, 45, 114,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3926, 3, 48, 120,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3944, 3, 51, 126,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3962, 3, 54, 132,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3980, 3, 57, 138,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3998, 3, 60, 144,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4016, 3, 63, 150,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4034, 3, 66, 156,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4052, 3, 84, 182,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4082, 3, 90, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4112, 3, 96, 202,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4142, 3, 102, 212,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4172, 3, 108, 222,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4202, 3, 114, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4232, 3, 120, 242,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4262, 3, 126, 252,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4292, 3, 132, 262,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4322, 3, 138, 272,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4352, 3, 144, 282,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4382, 3, 150, 292,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4412, 3, 182, 332,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4457, 3, 192, 347,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4502, 3, 202, 362,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4547, 3, 212, 377,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4592, 3, 222, 392,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4637, 3, 232, 407,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4682, 3, 242, 422,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4727, 3, 252, 437,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4772, 3, 262, 452,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4817, 3, 272, 467,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4862, 3, 282, 482,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4907, 3, 332, 539,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4970, 3, 347, 560,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5033, 3, 362, 581,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5096, 3, 377, 602,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5159, 3, 392, 623,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5222, 3, 407, 644,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5285, 3, 422, 665,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5348, 3, 437, 686,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5411, 3, 452, 707,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5474, 3, 467, 728,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5537, 3, 539, 805,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5621, 3, 560, 833,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5705, 3, 581, 861,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5789, 3, 602, 889,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5873, 3, 623, 917,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5957, 3, 644, 945,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6041, 3, 665, 973,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6125, 3, 686,
                                                                       1001, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6209, 3, 707,
                                                                       1029, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6293, 3, 805,
                                                                       1129, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6401, 3, 833,
                                                                       1165, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6509, 3, 861,
                                                                       1201, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6617, 3, 889,
                                                                       1237, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6725, 3, 917,
                                                                       1273, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6833, 3, 945,
                                                                       1309, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6941, 3, 973,
                                                                       1345, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7049, 3, 1001,
                                                                       1381, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7157, 3, 1129,
                                                                       1507, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7292, 3, 1165,
                                                                       1552, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7427, 3, 1201,
                                                                       1597, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7562, 3, 1237,
                                                                       1642, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7697, 3, 1273,
                                                                       1687, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7832, 3, 1309,
                                                                       1732, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7967, 3, 1345,
                                                                       1777, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8102, 3, 1507,
                                                                       1932, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8267, 3, 1552,
                                                                       1987, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8432, 3, 1597,
                                                                       2042, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8597, 3, 1642,
                                                                       2097, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8762, 3, 1687,
                                                                       2152, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8927, 3, 1732,
                                                                       2207, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9092, 3, 1932,
                                                                       2394, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9290, 3, 1987,
                                                                       2460, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9488, 3, 2042,
                                                                       2526, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9686, 3, 2097,
                                                                       2592, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9884, 3, 2152,
                                                                       2658, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 10082, 3, 2394,
                                                                       2880, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 10316, 3, 2460,
                                                                       2958, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 10550, 3, 2526,
                                                                       3036, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 10784, 3, 2592,
                                                                       3114, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 11018, 3, 2880,
                                                                       3374, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 11291, 3, 2958,
                                                                       3465, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 11564, 3, 3036,
                                                                       3556, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11837, 3, 7, 8,
                                                                       3647, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11843, 3, 8, 9,
                                                                       3650, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11849, 3, 9, 10,
                                                                       3653, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11855, 3, 10, 11,
                                                                       3656, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11861, 3, 11, 12,
                                                                       3659, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11867, 3, 12, 13,
                                                                       3662, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11873, 3, 13, 14,
                                                                       3665, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11879, 3, 14, 15,
                                                                       3668, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11885, 3, 15, 16,
                                                                       3671, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11891, 3, 16, 17,
                                                                       3674, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11897, 3, 17, 18,
                                                                       3677, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11903, 3, 18, 19,
                                                                       3680, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11909, 3, 19, 20,
                                                                       3683, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11915, 3, 20, 21,
                                                                       3686, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11921, 3, 21, 22,
                                                                       3689, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11927, 0, 3,
                                                                       11837, 3647, 11843, 3692,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11945, 0, 3,
                                                                       11843, 3650, 11849, 3701,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11963, 0, 3,
                                                                       11849, 3653, 11855, 3710,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11981, 0, 3,
                                                                       11855, 3656, 11861, 3719,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11999, 0, 3,
                                                                       11861, 3659, 11867, 3728,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12017, 0, 3,
                                                                       11867, 3662, 11873, 3737,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12035, 0, 3,
                                                                       11873, 3665, 11879, 3746,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12053, 0, 3,
                                                                       11879, 3668, 11885, 3755,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12071, 0, 3,
                                                                       11885, 3671, 11891, 3764,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12089, 0, 3,
                                                                       11891, 3674, 11897, 3773,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12107, 0, 3,
                                                                       11897, 3677, 11903, 3782,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12125, 0, 3,
                                                                       11903, 3680, 11909, 3791,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12143, 0, 3,
                                                                       11909, 3683, 11915, 3800,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12161, 0, 3,
                                                                       11915, 3686, 11921, 3809,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12179, 0, 3,
                                                                       11927, 3692, 11945, 72,
                                                                       78, 3818, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12215, 0, 3,
                                                                       11945, 3701, 11963, 78,
                                                                       84, 3836, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12251, 0, 3,
                                                                       11963, 3710, 11981, 84,
                                                                       90, 3854, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12287, 0, 3,
                                                                       11981, 3719, 11999, 90,
                                                                       96, 3872, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12323, 0, 3,
                                                                       11999, 3728, 12017, 96,
                                                                       102, 3890, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12359, 0, 3,
                                                                       12017, 3737, 12035, 102,
                                                                       108, 3908, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12395, 0, 3,
                                                                       12035, 3746, 12053, 108,
                                                                       114, 3926, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12431, 0, 3,
                                                                       12053, 3755, 12071, 114,
                                                                       120, 3944, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12467, 0, 3,
                                                                       12071, 3764, 12089, 120,
                                                                       126, 3962, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12503, 0, 3,
                                                                       12089, 3773, 12107, 126,
                                                                       132, 3980, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12539, 0, 3,
                                                                       12107, 3782, 12125, 132,
                                                                       138, 3998, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12575, 0, 3,
                                                                       12125, 3791, 12143, 138,
                                                                       144, 4016, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12611, 0, 3,
                                                                       12143, 3800, 12161, 144,
                                                                       150, 4034, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12647, 0, 3,
                                                                       12179, 3818, 12215, 162,
                                                                       172, 4052, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12707, 0, 3,
                                                                       12215, 3836, 12251, 172,
                                                                       182, 4082, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12767, 0, 3,
                                                                       12251, 3854, 12287, 182,
                                                                       192, 4112, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12827, 0, 3,
                                                                       12287, 3872, 12323, 192,
                                                                       202, 4142, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12887, 0, 3,
                                                                       12323, 3890, 12359, 202,
                                                                       212, 4172, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12947, 0, 3,
                                                                       12359, 3908, 12395, 212,
                                                                       222, 4202, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13007, 0, 3,
                                                                       12395, 3926, 12431, 222,
                                                                       232, 4232, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13067, 0, 3,
                                                                       12431, 3944, 12467, 232,
                                                                       242, 4262, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13127, 0, 3,
                                                                       12467, 3962, 12503, 242,
                                                                       252, 4292, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13187, 0, 3,
                                                                       12503, 3980, 12539, 252,
                                                                       262, 4322, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13247, 0, 3,
                                                                       12539, 3998, 12575, 262,
                                                                       272, 4352, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13307, 0, 3,
                                                                       12575, 4016, 12611, 272,
                                                                       282, 4382, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13367, 0, 3,
                                                                       12647, 4052, 12707, 302,
                                                                       317, 4412, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13457, 0, 3,
                                                                       12707, 4082, 12767, 317,
                                                                       332, 4457, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13547, 0, 3,
                                                                       12767, 4112, 12827, 332,
                                                                       347, 4502, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13637, 0, 3,
                                                                       12827, 4142, 12887, 347,
                                                                       362, 4547, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13727, 0, 3,
                                                                       12887, 4172, 12947, 362,
                                                                       377, 4592, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13817, 0, 3,
                                                                       12947, 4202, 13007, 377,
                                                                       392, 4637, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13907, 0, 3,
                                                                       13007, 4232, 13067, 392,
                                                                       407, 4682, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13997, 0, 3,
                                                                       13067, 4262, 13127, 407,
                                                                       422, 4727, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14087, 0, 3,
                                                                       13127, 4292, 13187, 422,
                                                                       437, 4772, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14177, 0, 3,
                                                                       13187, 4322, 13247, 437,
                                                                       452, 4817, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14267, 0, 3,
                                                                       13247, 4352, 13307, 452,
                                                                       467, 4862, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14357, 0, 3,
                                                                       13367, 4412, 13457, 497,
                                                                       518, 4907, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14483, 0, 3,
                                                                       13457, 4457, 13547, 518,
                                                                       539, 4970, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14609, 0, 3,
                                                                       13547, 4502, 13637, 539,
                                                                       560, 5033, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14735, 0, 3,
                                                                       13637, 4547, 13727, 560,
                                                                       581, 5096, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14861, 0, 3,
                                                                       13727, 4592, 13817, 581,
                                                                       602, 5159, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14987, 0, 3,
                                                                       13817, 4637, 13907, 602,
                                                                       623, 5222, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15113, 0, 3,
                                                                       13907, 4682, 13997, 623,
                                                                       644, 5285, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15239, 0, 3,
                                                                       13997, 4727, 14087, 644,
                                                                       665, 5348, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15365, 0, 3,
                                                                       14087, 4772, 14177, 665,
                                                                       686, 5411, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15491, 0, 3,
                                                                       14177, 4817, 14267, 686,
                                                                       707, 5474, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15617, 0, 3,
                                                                       14357, 4907, 14483, 749,
                                                                       777, 5537, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15785, 0, 3,
                                                                       14483, 4970, 14609, 777,
                                                                       805, 5621, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15953, 0, 3,
                                                                       14609, 5033, 14735, 805,
                                                                       833, 5705, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16121, 0, 3,
                                                                       14735, 5096, 14861, 833,
                                                                       861, 5789, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16289, 0, 3,
                                                                       14861, 5159, 14987, 861,
                                                                       889, 5873, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16457, 0, 3,
                                                                       14987, 5222, 15113, 889,
                                                                       917, 5957, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16625, 0, 3,
                                                                       15113, 5285, 15239, 917,
                                                                       945, 6041, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16793, 0, 3,
                                                                       15239, 5348, 15365, 945,
                                                                       973, 6125, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16961, 0, 3,
                                                                       15365, 5411, 15491, 973,
                                                                       1001, 6209, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17129, 0, 3,
                                                                       15617, 5537, 15785, 1057,
                                                                       1093, 6293, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17345, 0, 3,
                                                                       15785, 5621, 15953, 1093,
                                                                       1129, 6401, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17561, 0, 3,
                                                                       15953, 5705, 16121, 1129,
                                                                       1165, 6509, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17777, 0, 3,
                                                                       16121, 5789, 16289, 1165,
                                                                       1201, 6617, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17993, 0, 3,
                                                                       16289, 5873, 16457, 1201,
                                                                       1237, 6725, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18209, 0, 3,
                                                                       16457, 5957, 16625, 1237,
                                                                       1273, 6833, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18425, 0, 3,
                                                                       16625, 6041, 16793, 1273,
                                                                       1309, 6941, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18641, 0, 3,
                                                                       16793, 6125, 16961, 1309,
                                                                       1345, 7049, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18857, 0, 3,
                                                                       17129, 6293, 17345, 1417,
                                                                       1462, 7157, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 19127, 0, 3,
                                                                       17345, 6401, 17561, 1462,
                                                                       1507, 7292, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 19397, 0, 3,
                                                                       17561, 6509, 17777, 1507,
                                                                       1552, 7427, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 19667, 0, 3,
                                                                       17777, 6617, 17993, 1552,
                                                                       1597, 7562, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 19937, 0, 3,
                                                                       17993, 6725, 18209, 1597,
                                                                       1642, 7697, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 20207, 0, 3,
                                                                       18209, 6833, 18425, 1642,
                                                                       1687, 7832, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 20477, 0, 3,
                                                                       18425, 6941, 18641, 1687,
                                                                       1732, 7967, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 20747, 0, 3,
                                                                       18857, 7157, 19127, 1822,
                                                                       1877, 8102, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 21077, 0, 3,
                                                                       19127, 7292, 19397, 1877,
                                                                       1932, 8267, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 21407, 0, 3,
                                                                       19397, 7427, 19667, 1932,
                                                                       1987, 8432, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 21737, 0, 3,
                                                                       19667, 7562, 19937, 1987,
                                                                       2042, 8597, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 22067, 0, 3,
                                                                       19937, 7697, 20207, 2042,
                                                                       2097, 8762, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 22397, 0, 3,
                                                                       20207, 7832, 20477, 2097,
                                                                       2152, 8927, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 22727, 0, 3,
                                                                       20747, 8102, 21077, 2262,
                                                                       2328, 9092, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 23123, 0, 3,
                                                                       21077, 8267, 21407, 2328,
                                                                       2394, 9290, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 23519, 0, 3,
                                                                       21407, 8432, 21737, 2394,
                                                                       2460, 9488, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 23915, 0, 3,
                                                                       21737, 8597, 22067, 2460,
                                                                       2526, 9686, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 24311, 0, 3,
                                                                       22067, 8762, 22397, 2526,
                                                                       2592, 9884, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 24707, 0, 3,
                                                                       22727, 9092, 23123, 2724,
                                                                       2802, 10082, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 25175, 0, 3,
                                                                       23123, 9290, 23519, 2802,
                                                                       2880, 10316, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 25643, 0, 3,
                                                                       23519, 9488, 23915, 2880,
                                                                       2958, 10550, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 26111, 0, 3,
                                                                       23915, 9686, 24311, 2958,
                                                                       3036, 10784, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 26579, 0, 3,
                                                                       24707, 10082, 25175, 3192,
                                                                       3283, 11018, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 27125, 0, 3,
                                                                       25175, 10316, 25643, 3283,
                                                                       3374, 11291, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 27671, 0, 3,
                                                                       25643, 10550, 26111, 3374,
                                                                       3465, 11564, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28217, 3, 3647,
                                                                       3650, 11849, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28227, 3, 3650,
                                                                       3653, 11855, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28237, 3, 3653,
                                                                       3656, 11861, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28247, 3, 3656,
                                                                       3659, 11867, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28257, 3, 3659,
                                                                       3662, 11873, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28267, 3, 3662,
                                                                       3665, 11879, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28277, 3, 3665,
                                                                       3668, 11885, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28287, 3, 3668,
                                                                       3671, 11891, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28297, 3, 3671,
                                                                       3674, 11897, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28307, 3, 3674,
                                                                       3677, 11903, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28317, 3, 3677,
                                                                       3680, 11909, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28327, 3, 3680,
                                                                       3683, 11915, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28337, 3, 3683,
                                                                       3686, 11921, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28347, 0, 3,
                                                                       28217, 11849, 28227, 3692,
                                                                       3701, 11963, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28377, 0, 3,
                                                                       28227, 11855, 28237, 3701,
                                                                       3710, 11981, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28407, 0, 3,
                                                                       28237, 11861, 28247, 3710,
                                                                       3719, 11999, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28437, 0, 3,
                                                                       28247, 11867, 28257, 3719,
                                                                       3728, 12017, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28467, 0, 3,
                                                                       28257, 11873, 28267, 3728,
                                                                       3737, 12035, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28497, 0, 3,
                                                                       28267, 11879, 28277, 3737,
                                                                       3746, 12053, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28527, 0, 3,
                                                                       28277, 11885, 28287, 3746,
                                                                       3755, 12071, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28557, 0, 3,
                                                                       28287, 11891, 28297, 3755,
                                                                       3764, 12089, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28587, 0, 3,
                                                                       28297, 11897, 28307, 3764,
                                                                       3773, 12107, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28617, 0, 3,
                                                                       28307, 11903, 28317, 3773,
                                                                       3782, 12125, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28647, 0, 3,
                                                                       28317, 11909, 28327, 3782,
                                                                       3791, 12143, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28677, 0, 3,
                                                                       28327, 11915, 28337, 3791,
                                                                       3800, 12161, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28707, 0, 3,
                                                                       28347, 11963, 28377, 3818,
                                                                       3836, 12251, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28767, 0, 3,
                                                                       28377, 11981, 28407, 3836,
                                                                       3854, 12287, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28827, 0, 3,
                                                                       28407, 11999, 28437, 3854,
                                                                       3872, 12323, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28887, 0, 3,
                                                                       28437, 12017, 28467, 3872,
                                                                       3890, 12359, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28947, 0, 3,
                                                                       28467, 12035, 28497, 3890,
                                                                       3908, 12395, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29007, 0, 3,
                                                                       28497, 12053, 28527, 3908,
                                                                       3926, 12431, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29067, 0, 3,
                                                                       28527, 12071, 28557, 3926,
                                                                       3944, 12467, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29127, 0, 3,
                                                                       28557, 12089, 28587, 3944,
                                                                       3962, 12503, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29187, 0, 3,
                                                                       28587, 12107, 28617, 3962,
                                                                       3980, 12539, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29247, 0, 3,
                                                                       28617, 12125, 28647, 3980,
                                                                       3998, 12575, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29307, 0, 3,
                                                                       28647, 12143, 28677, 3998,
                                                                       4016, 12611, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29367, 0, 3,
                                                                       28707, 12251, 28767, 4052,
                                                                       4082, 12767, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29467, 0, 3,
                                                                       28767, 12287, 28827, 4082,
                                                                       4112, 12827, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29567, 0, 3,
                                                                       28827, 12323, 28887, 4112,
                                                                       4142, 12887, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29667, 0, 3,
                                                                       28887, 12359, 28947, 4142,
                                                                       4172, 12947, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29767, 0, 3,
                                                                       28947, 12395, 29007, 4172,
                                                                       4202, 13007, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29867, 0, 3,
                                                                       29007, 12431, 29067, 4202,
                                                                       4232, 13067, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29967, 0, 3,
                                                                       29067, 12467, 29127, 4232,
                                                                       4262, 13127, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30067, 0, 3,
                                                                       29127, 12503, 29187, 4262,
                                                                       4292, 13187, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30167, 0, 3,
                                                                       29187, 12539, 29247, 4292,
                                                                       4322, 13247, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30267, 0, 3,
                                                                       29247, 12575, 29307, 4322,
                                                                       4352, 13307, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30367, 0, 3,
                                                                       29367, 12767, 29467, 4412,
                                                                       4457, 13547, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30517, 0, 3,
                                                                       29467, 12827, 29567, 4457,
                                                                       4502, 13637, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30667, 0, 3,
                                                                       29567, 12887, 29667, 4502,
                                                                       4547, 13727, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30817, 0, 3,
                                                                       29667, 12947, 29767, 4547,
                                                                       4592, 13817, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30967, 0, 3,
                                                                       29767, 13007, 29867, 4592,
                                                                       4637, 13907, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31117, 0, 3,
                                                                       29867, 13067, 29967, 4637,
                                                                       4682, 13997, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31267, 0, 3,
                                                                       29967, 13127, 30067, 4682,
                                                                       4727, 14087, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31417, 0, 3,
                                                                       30067, 13187, 30167, 4727,
                                                                       4772, 14177, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31567, 0, 3,
                                                                       30167, 13247, 30267, 4772,
                                                                       4817, 14267, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 31717, 0, 3,
                                                                       30367, 13547, 30517, 4907,
                                                                       4970, 14609, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 31927, 0, 3,
                                                                       30517, 13637, 30667, 4970,
                                                                       5033, 14735, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32137, 0, 3,
                                                                       30667, 13727, 30817, 5033,
                                                                       5096, 14861, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32347, 0, 3,
                                                                       30817, 13817, 30967, 5096,
                                                                       5159, 14987, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32557, 0, 3,
                                                                       30967, 13907, 31117, 5159,
                                                                       5222, 15113, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32767, 0, 3,
                                                                       31117, 13997, 31267, 5222,
                                                                       5285, 15239, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32977, 0, 3,
                                                                       31267, 14087, 31417, 5285,
                                                                       5348, 15365, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33187, 0, 3,
                                                                       31417, 14177, 31567, 5348,
                                                                       5411, 15491, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 33397, 0, 3,
                                                                       31717, 14609, 31927, 5537,
                                                                       5621, 15953, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 33677, 0, 3,
                                                                       31927, 14735, 32137, 5621,
                                                                       5705, 16121, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 33957, 0, 3,
                                                                       32137, 14861, 32347, 5705,
                                                                       5789, 16289, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 34237, 0, 3,
                                                                       32347, 14987, 32557, 5789,
                                                                       5873, 16457, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 34517, 0, 3,
                                                                       32557, 15113, 32767, 5873,
                                                                       5957, 16625, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 34797, 0, 3,
                                                                       32767, 15239, 32977, 5957,
                                                                       6041, 16793, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 35077, 0, 3,
                                                                       32977, 15365, 33187, 6041,
                                                                       6125, 16961, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 35357, 0, 3,
                                                                       33397, 15953, 33677, 6293,
                                                                       6401, 17561, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 35717, 0, 3,
                                                                       33677, 16121, 33957, 6401,
                                                                       6509, 17777, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 36077, 0, 3,
                                                                       33957, 16289, 34237, 6509,
                                                                       6617, 17993, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 36437, 0, 3,
                                                                       34237, 16457, 34517, 6617,
                                                                       6725, 18209, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 36797, 0, 3,
                                                                       34517, 16625, 34797, 6725,
                                                                       6833, 18425, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 37157, 0, 3,
                                                                       34797, 16793, 35077, 6833,
                                                                       6941, 18641, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 37517, 0, 3,
                                                                       35357, 17561, 35717, 7157,
                                                                       7292, 19397, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 37967, 0, 3,
                                                                       35717, 17777, 36077, 7292,
                                                                       7427, 19667, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 38417, 0, 3,
                                                                       36077, 17993, 36437, 7427,
                                                                       7562, 19937, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 38867, 0, 3,
                                                                       36437, 18209, 36797, 7562,
                                                                       7697, 20207, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 39317, 0, 3,
                                                                       36797, 18425, 37157, 7697,
                                                                       7832, 20477, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 39767, 0, 3,
                                                                       37517, 19397, 37967, 8102,
                                                                       8267, 21407, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 40317, 0, 3,
                                                                       37967, 19667, 38417, 8267,
                                                                       8432, 21737, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 40867, 0, 3,
                                                                       38417, 19937, 38867, 8432,
                                                                       8597, 22067, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 41417, 0, 3,
                                                                       38867, 20207, 39317, 8597,
                                                                       8762, 22397, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 41967, 0, 3,
                                                                       39767, 21407, 40317, 9092,
                                                                       9290, 23519, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 42627, 0, 3,
                                                                       40317, 21737, 40867, 9290,
                                                                       9488, 23915, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 43287, 0, 3,
                                                                       40867, 22067, 41417, 9488,
                                                                       9686, 24311, ncols, gamma,
                                                                       p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 43947, 0, 3,
                                                                       41967, 23519, 42627,
                                                                       10082, 10316, 25643,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 44727, 0, 3,
                                                                       42627, 23915, 43287,
                                                                       10316, 10550, 26111,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 45507, 0, 3,
                                                                       43947, 25643, 44727,
                                                                       11018, 11291, 27671,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46417, 3, 11837,
                                                                       11843, 28217, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46432, 3, 11843,
                                                                       11849, 28227, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46447, 3, 11849,
                                                                       11855, 28237, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46462, 3, 11855,
                                                                       11861, 28247, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46477, 3, 11861,
                                                                       11867, 28257, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46492, 3, 11867,
                                                                       11873, 28267, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46507, 3, 11873,
                                                                       11879, 28277, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46522, 3, 11879,
                                                                       11885, 28287, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46537, 3, 11885,
                                                                       11891, 28297, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46552, 3, 11891,
                                                                       11897, 28307, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46567, 3, 11897,
                                                                       11903, 28317, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46582, 3, 11903,
                                                                       11909, 28327, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 46597, 3, 11909,
                                                                       11915, 28337, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 46612, 0, 3,
                                                                       46417, 28217, 46432,
                                                                       11927, 11945, 28347,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 46657, 0, 3,
                                                                       46432, 28227, 46447,
                                                                       11945, 11963, 28377,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 46702, 0, 3,
                                                                       46447, 28237, 46462,
                                                                       11963, 11981, 28407,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 46747, 0, 3,
                                                                       46462, 28247, 46477,
                                                                       11981, 11999, 28437,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 46792, 0, 3,
                                                                       46477, 28257, 46492,
                                                                       11999, 12017, 28467,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 46837, 0, 3,
                                                                       46492, 28267, 46507,
                                                                       12017, 12035, 28497,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 46882, 0, 3,
                                                                       46507, 28277, 46522,
                                                                       12035, 12053, 28527,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 46927, 0, 3,
                                                                       46522, 28287, 46537,
                                                                       12053, 12071, 28557,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 46972, 0, 3,
                                                                       46537, 28297, 46552,
                                                                       12071, 12089, 28587,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 47017, 0, 3,
                                                                       46552, 28307, 46567,
                                                                       12089, 12107, 28617,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 47062, 0, 3,
                                                                       46567, 28317, 46582,
                                                                       12107, 12125, 28647,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 47107, 0, 3,
                                                                       46582, 28327, 46597,
                                                                       12125, 12143, 28677,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 47152, 0, 3,
                                                                       46612, 28347, 46657,
                                                                       12179, 12215, 28707,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 47242, 0, 3,
                                                                       46657, 28377, 46702,
                                                                       12215, 12251, 28767,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 47332, 0, 3,
                                                                       46702, 28407, 46747,
                                                                       12251, 12287, 28827,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 47422, 0, 3,
                                                                       46747, 28437, 46792,
                                                                       12287, 12323, 28887,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 47512, 0, 3,
                                                                       46792, 28467, 46837,
                                                                       12323, 12359, 28947,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 47602, 0, 3,
                                                                       46837, 28497, 46882,
                                                                       12359, 12395, 29007,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 47692, 0, 3,
                                                                       46882, 28527, 46927,
                                                                       12395, 12431, 29067,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 47782, 0, 3,
                                                                       46927, 28557, 46972,
                                                                       12431, 12467, 29127,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 47872, 0, 3,
                                                                       46972, 28587, 47017,
                                                                       12467, 12503, 29187,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 47962, 0, 3,
                                                                       47017, 28617, 47062,
                                                                       12503, 12539, 29247,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 48052, 0, 3,
                                                                       47062, 28647, 47107,
                                                                       12539, 12575, 29307,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 48142, 0, 3,
                                                                       47152, 28707, 47242,
                                                                       12647, 12707, 29367,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 48292, 0, 3,
                                                                       47242, 28767, 47332,
                                                                       12707, 12767, 29467,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 48442, 0, 3,
                                                                       47332, 28827, 47422,
                                                                       12767, 12827, 29567,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 48592, 0, 3,
                                                                       47422, 28887, 47512,
                                                                       12827, 12887, 29667,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 48742, 0, 3,
                                                                       47512, 28947, 47602,
                                                                       12887, 12947, 29767,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 48892, 0, 3,
                                                                       47602, 29007, 47692,
                                                                       12947, 13007, 29867,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 49042, 0, 3,
                                                                       47692, 29067, 47782,
                                                                       13007, 13067, 29967,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 49192, 0, 3,
                                                                       47782, 29127, 47872,
                                                                       13067, 13127, 30067,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 49342, 0, 3,
                                                                       47872, 29187, 47962,
                                                                       13127, 13187, 30167,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 49492, 0, 3,
                                                                       47962, 29247, 48052,
                                                                       13187, 13247, 30267,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 49642, 0, 3,
                                                                       48142, 29367, 48292,
                                                                       13367, 13457, 30367,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 49867, 0, 3,
                                                                       48292, 29467, 48442,
                                                                       13457, 13547, 30517,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 50092, 0, 3,
                                                                       48442, 29567, 48592,
                                                                       13547, 13637, 30667,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 50317, 0, 3,
                                                                       48592, 29667, 48742,
                                                                       13637, 13727, 30817,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 50542, 0, 3,
                                                                       48742, 29767, 48892,
                                                                       13727, 13817, 30967,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 50767, 0, 3,
                                                                       48892, 29867, 49042,
                                                                       13817, 13907, 31117,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 50992, 0, 3,
                                                                       49042, 29967, 49192,
                                                                       13907, 13997, 31267,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 51217, 0, 3,
                                                                       49192, 30067, 49342,
                                                                       13997, 14087, 31417,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 51442, 0, 3,
                                                                       49342, 30167, 49492,
                                                                       14087, 14177, 31567,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 51667, 0, 3,
                                                                       49642, 30367, 49867,
                                                                       14357, 14483, 31717,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 51982, 0, 3,
                                                                       49867, 30517, 50092,
                                                                       14483, 14609, 31927,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 52297, 0, 3,
                                                                       50092, 30667, 50317,
                                                                       14609, 14735, 32137,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 52612, 0, 3,
                                                                       50317, 30817, 50542,
                                                                       14735, 14861, 32347,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 52927, 0, 3,
                                                                       50542, 30967, 50767,
                                                                       14861, 14987, 32557,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 53242, 0, 3,
                                                                       50767, 31117, 50992,
                                                                       14987, 15113, 32767,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 53557, 0, 3,
                                                                       50992, 31267, 51217,
                                                                       15113, 15239, 32977,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 53872, 0, 3,
                                                                       51217, 31417, 51442,
                                                                       15239, 15365, 33187,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 54187, 0, 3,
                                                                       51667, 31717, 51982,
                                                                       15617, 15785, 33397,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 54607, 0, 3,
                                                                       51982, 31927, 52297,
                                                                       15785, 15953, 33677,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 55027, 0, 3,
                                                                       52297, 32137, 52612,
                                                                       15953, 16121, 33957,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 55447, 0, 3,
                                                                       52612, 32347, 52927,
                                                                       16121, 16289, 34237,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 55867, 0, 3,
                                                                       52927, 32557, 53242,
                                                                       16289, 16457, 34517,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 56287, 0, 3,
                                                                       53242, 32767, 53557,
                                                                       16457, 16625, 34797,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 56707, 0, 3,
                                                                       53557, 32977, 53872,
                                                                       16625, 16793, 35077,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 57127, 0, 3,
                                                                       54187, 33397, 54607,
                                                                       17129, 17345, 35357,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 57667, 0, 3,
                                                                       54607, 33677, 55027,
                                                                       17345, 17561, 35717,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 58207, 0, 3,
                                                                       55027, 33957, 55447,
                                                                       17561, 17777, 36077,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 58747, 0, 3,
                                                                       55447, 34237, 55867,
                                                                       17777, 17993, 36437,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 59287, 0, 3,
                                                                       55867, 34517, 56287,
                                                                       17993, 18209, 36797,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 59827, 0, 3,
                                                                       56287, 34797, 56707,
                                                                       18209, 18425, 37157,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 60367, 0, 3,
                                                                       57127, 35357, 57667,
                                                                       18857, 19127, 37517,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 61042, 0, 3,
                                                                       57667, 35717, 58207,
                                                                       19127, 19397, 37967,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 61717, 0, 3,
                                                                       58207, 36077, 58747,
                                                                       19397, 19667, 38417,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 62392, 0, 3,
                                                                       58747, 36437, 59287,
                                                                       19667, 19937, 38867,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 63067, 0, 3,
                                                                       59287, 36797, 59827,
                                                                       19937, 20207, 39317,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 63742, 0, 3,
                                                                       60367, 37517, 61042,
                                                                       20747, 21077, 39767,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 64567, 0, 3,
                                                                       61042, 37967, 61717,
                                                                       21077, 21407, 40317,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 65392, 0, 3,
                                                                       61717, 38417, 62392,
                                                                       21407, 21737, 40867,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 66217, 0, 3,
                                                                       62392, 38867, 63067,
                                                                       21737, 22067, 41417,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 67042, 0, 3,
                                                                       63742, 39767, 64567,
                                                                       22727, 23123, 41967,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 68032, 0, 3,
                                                                       64567, 40317, 65392,
                                                                       23123, 23519, 42627,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 69022, 0, 3,
                                                                       65392, 40867, 66217,
                                                                       23519, 23915, 43287,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 70012, 0, 3,
                                                                       67042, 41967, 68032,
                                                                       24707, 25175, 43947,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 71182, 0, 3,
                                                                       68032, 42627, 69022,
                                                                       25175, 25643, 44727,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 72352, 0, 3,
                                                                       70012, 43947, 71182,
                                                                       26579, 27125, 45507,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 73717, 54187, 420, ncols);

                    simdfunc::contract_primitives(buffer, 74389, 57127, 540, ncols);

                    simdfunc::contract_primitives(buffer, 75253, 60367, 675, ncols);

                    simdfunc::contract_primitives(buffer, 76333, 63742, 825, ncols);

                    simdfunc::contract_primitives(buffer, 77653, 67042, 990, ncols);

                    simdfunc::contract_primitives(buffer, 79237, 70012, 1170, ncols);

                    simdfunc::contract_primitives(buffer, 81109, 72352, 1365, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 74137, 73717, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 74929, 74389, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 75928, 75253, 45, 1, nmax);

        simdtrf::transform_g_inner(buffer, 77158, 76333, 55, 1, nmax);

        simdtrf::transform_g_inner(buffer, 78643, 77653, 66, 1, nmax);

        simdtrf::transform_g_inner(buffer, 80407, 79237, 78, 1, nmax);

        simdtrf::transform_g_inner(buffer, 82474, 81109, 91, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 83293, 74137, 74929, 9, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 84049, 74929, 75928, 9, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 85021, 75928, 77158, 9, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 86236, 77158, 78643, 9, nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 87721, 78643, 80407, 9, nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 89503, 80407, 82474, 9, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 91609, 83293, 84049, 9, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 93121, 84049, 85021, 9, nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 95065, 85021, 86236, 9, nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 97495, 86236, 87721, 9, nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 100465, 87721, 89503, 9,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 104029, 91609, 93121, 9,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 106549, 93121, 95065, 9,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 109789, 95065, 97495, 9,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 113839, 97495, 100465, 9,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 118789, 104029, 106549, 9,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 122569, 106549, 109789, 9,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 127429, 109789, 113839, 9,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 133504, 118789, 122569, 9,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 138796, 122569, 127429, 9,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 145600, 133504, 138796, 9,
                                             nmax);

        simdtrf::transform_i_inner(buffer, 152656, 145600, 28, 9, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 152656, 117, nmax);
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
