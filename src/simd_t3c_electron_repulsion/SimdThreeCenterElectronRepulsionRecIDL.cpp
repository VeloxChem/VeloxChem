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


#include "SimdThreeCenterElectronRepulsionRecIDL.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_idl_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_idl_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 139840, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1105 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 139840, 124582, 5993, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1822, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1825, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1828, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1831, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1834, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1837, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1840, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1843, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1846, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1849, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1852, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1855, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1858, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1861, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1864, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1867, 3, 9, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1876, 3, 10, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1885, 3, 11, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1894, 3, 12, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1903, 3, 13, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1912, 3, 14, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1921, 3, 15, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1930, 3, 16, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1939, 3, 17, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1948, 3, 18, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1957, 3, 19, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1966, 3, 20, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1975, 3, 21, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1984, 3, 22, 69,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1993, 3, 30, 84,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2011, 3, 33, 90,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2029, 3, 36, 96,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2047, 3, 39, 102,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2065, 3, 42, 108,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2083, 3, 45, 114,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2101, 3, 48, 120,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2119, 3, 51, 126,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2137, 3, 54, 132,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2155, 3, 57, 138,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2173, 3, 60, 144,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2191, 3, 63, 150,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2209, 3, 66, 156,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2227, 3, 84, 182,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2257, 3, 90, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2287, 3, 96, 202,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2317, 3, 102, 212,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2347, 3, 108, 222,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2377, 3, 114, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2407, 3, 120, 242,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2437, 3, 126, 252,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2467, 3, 132, 262,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2497, 3, 138, 272,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2527, 3, 144, 282,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2557, 3, 150, 292,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2587, 3, 182, 332,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2632, 3, 192, 347,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2677, 3, 202, 362,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2722, 3, 212, 377,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2767, 3, 222, 392,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2812, 3, 232, 407,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2857, 3, 242, 422,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2902, 3, 252, 437,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2947, 3, 262, 452,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2992, 3, 272, 467,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3037, 3, 282, 482,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3082, 3, 332, 539,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3145, 3, 347, 560,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3208, 3, 362, 581,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3271, 3, 377, 602,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3334, 3, 392, 623,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3397, 3, 407, 644,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3460, 3, 422, 665,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3523, 3, 437, 686,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3586, 3, 452, 707,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3649, 3, 467, 728,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3712, 3, 539, 805,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3796, 3, 560, 833,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3880, 3, 581, 861,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3964, 3, 602, 889,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4048, 3, 623, 917,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4132, 3, 644, 945,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4216, 3, 665, 973,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4300, 3, 686,
                                                                       1001, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4384, 3, 707,
                                                                       1029, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4468, 3, 805,
                                                                       1129, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4576, 3, 833,
                                                                       1165, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4684, 3, 861,
                                                                       1201, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4792, 3, 889,
                                                                       1237, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4900, 3, 917,
                                                                       1273, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5008, 3, 945,
                                                                       1309, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5116, 3, 973,
                                                                       1345, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5224, 3, 1001,
                                                                       1381, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5332, 3, 1129,
                                                                       1507, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5467, 3, 1165,
                                                                       1552, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5602, 3, 1201,
                                                                       1597, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5737, 3, 1237,
                                                                       1642, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5872, 3, 1273,
                                                                       1687, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6007, 3, 1309,
                                                                       1732, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6142, 3, 1345,
                                                                       1777, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6277, 3, 7, 8,
                                                                       1822, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6283, 3, 8, 9,
                                                                       1825, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6289, 3, 9, 10,
                                                                       1828, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6295, 3, 10, 11,
                                                                       1831, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6301, 3, 11, 12,
                                                                       1834, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6307, 3, 12, 13,
                                                                       1837, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6313, 3, 13, 14,
                                                                       1840, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6319, 3, 14, 15,
                                                                       1843, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6325, 3, 15, 16,
                                                                       1846, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6331, 3, 16, 17,
                                                                       1849, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6337, 3, 17, 18,
                                                                       1852, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6343, 3, 18, 19,
                                                                       1855, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6349, 3, 19, 20,
                                                                       1858, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6355, 3, 20, 21,
                                                                       1861, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6361, 3, 21, 22,
                                                                       1864, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6367, 0, 3, 6277,
                                                                       1822, 6283, 1867, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6385, 0, 3, 6283,
                                                                       1825, 6289, 1876, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6403, 0, 3, 6289,
                                                                       1828, 6295, 1885, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6421, 0, 3, 6295,
                                                                       1831, 6301, 1894, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6439, 0, 3, 6301,
                                                                       1834, 6307, 1903, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6457, 0, 3, 6307,
                                                                       1837, 6313, 1912, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6475, 0, 3, 6313,
                                                                       1840, 6319, 1921, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6493, 0, 3, 6319,
                                                                       1843, 6325, 1930, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6511, 0, 3, 6325,
                                                                       1846, 6331, 1939, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6529, 0, 3, 6331,
                                                                       1849, 6337, 1948, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6547, 0, 3, 6337,
                                                                       1852, 6343, 1957, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6565, 0, 3, 6343,
                                                                       1855, 6349, 1966, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6583, 0, 3, 6349,
                                                                       1858, 6355, 1975, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6601, 0, 3, 6355,
                                                                       1861, 6361, 1984, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6619, 0, 3, 6367,
                                                                       1867, 6385, 72, 78, 1993,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6655, 0, 3, 6385,
                                                                       1876, 6403, 78, 84, 2011,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6691, 0, 3, 6403,
                                                                       1885, 6421, 84, 90, 2029,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6727, 0, 3, 6421,
                                                                       1894, 6439, 90, 96, 2047,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6763, 0, 3, 6439,
                                                                       1903, 6457, 96, 102, 2065,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6799, 0, 3, 6457,
                                                                       1912, 6475, 102, 108,
                                                                       2083, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6835, 0, 3, 6475,
                                                                       1921, 6493, 108, 114,
                                                                       2101, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6871, 0, 3, 6493,
                                                                       1930, 6511, 114, 120,
                                                                       2119, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6907, 0, 3, 6511,
                                                                       1939, 6529, 120, 126,
                                                                       2137, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6943, 0, 3, 6529,
                                                                       1948, 6547, 126, 132,
                                                                       2155, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6979, 0, 3, 6547,
                                                                       1957, 6565, 132, 138,
                                                                       2173, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7015, 0, 3, 6565,
                                                                       1966, 6583, 138, 144,
                                                                       2191, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7051, 0, 3, 6583,
                                                                       1975, 6601, 144, 150,
                                                                       2209, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7087, 0, 3, 6619,
                                                                       1993, 6655, 162, 172,
                                                                       2227, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7147, 0, 3, 6655,
                                                                       2011, 6691, 172, 182,
                                                                       2257, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7207, 0, 3, 6691,
                                                                       2029, 6727, 182, 192,
                                                                       2287, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7267, 0, 3, 6727,
                                                                       2047, 6763, 192, 202,
                                                                       2317, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7327, 0, 3, 6763,
                                                                       2065, 6799, 202, 212,
                                                                       2347, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7387, 0, 3, 6799,
                                                                       2083, 6835, 212, 222,
                                                                       2377, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7447, 0, 3, 6835,
                                                                       2101, 6871, 222, 232,
                                                                       2407, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7507, 0, 3, 6871,
                                                                       2119, 6907, 232, 242,
                                                                       2437, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7567, 0, 3, 6907,
                                                                       2137, 6943, 242, 252,
                                                                       2467, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7627, 0, 3, 6943,
                                                                       2155, 6979, 252, 262,
                                                                       2497, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7687, 0, 3, 6979,
                                                                       2173, 7015, 262, 272,
                                                                       2527, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7747, 0, 3, 7015,
                                                                       2191, 7051, 272, 282,
                                                                       2557, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7807, 0, 3, 7087,
                                                                       2227, 7147, 302, 317,
                                                                       2587, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7897, 0, 3, 7147,
                                                                       2257, 7207, 317, 332,
                                                                       2632, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7987, 0, 3, 7207,
                                                                       2287, 7267, 332, 347,
                                                                       2677, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8077, 0, 3, 7267,
                                                                       2317, 7327, 347, 362,
                                                                       2722, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8167, 0, 3, 7327,
                                                                       2347, 7387, 362, 377,
                                                                       2767, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8257, 0, 3, 7387,
                                                                       2377, 7447, 377, 392,
                                                                       2812, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8347, 0, 3, 7447,
                                                                       2407, 7507, 392, 407,
                                                                       2857, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8437, 0, 3, 7507,
                                                                       2437, 7567, 407, 422,
                                                                       2902, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8527, 0, 3, 7567,
                                                                       2467, 7627, 422, 437,
                                                                       2947, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8617, 0, 3, 7627,
                                                                       2497, 7687, 437, 452,
                                                                       2992, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8707, 0, 3, 7687,
                                                                       2527, 7747, 452, 467,
                                                                       3037, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8797, 0, 3, 7807,
                                                                       2587, 7897, 497, 518,
                                                                       3082, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8923, 0, 3, 7897,
                                                                       2632, 7987, 518, 539,
                                                                       3145, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9049, 0, 3, 7987,
                                                                       2677, 8077, 539, 560,
                                                                       3208, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9175, 0, 3, 8077,
                                                                       2722, 8167, 560, 581,
                                                                       3271, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9301, 0, 3, 8167,
                                                                       2767, 8257, 581, 602,
                                                                       3334, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9427, 0, 3, 8257,
                                                                       2812, 8347, 602, 623,
                                                                       3397, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9553, 0, 3, 8347,
                                                                       2857, 8437, 623, 644,
                                                                       3460, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9679, 0, 3, 8437,
                                                                       2902, 8527, 644, 665,
                                                                       3523, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9805, 0, 3, 8527,
                                                                       2947, 8617, 665, 686,
                                                                       3586, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9931, 0, 3, 8617,
                                                                       2992, 8707, 686, 707,
                                                                       3649, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10057, 0, 3, 8797,
                                                                       3082, 8923, 749, 777,
                                                                       3712, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10225, 0, 3, 8923,
                                                                       3145, 9049, 777, 805,
                                                                       3796, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10393, 0, 3, 9049,
                                                                       3208, 9175, 805, 833,
                                                                       3880, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10561, 0, 3, 9175,
                                                                       3271, 9301, 833, 861,
                                                                       3964, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10729, 0, 3, 9301,
                                                                       3334, 9427, 861, 889,
                                                                       4048, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10897, 0, 3, 9427,
                                                                       3397, 9553, 889, 917,
                                                                       4132, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11065, 0, 3, 9553,
                                                                       3460, 9679, 917, 945,
                                                                       4216, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11233, 0, 3, 9679,
                                                                       3523, 9805, 945, 973,
                                                                       4300, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11401, 0, 3, 9805,
                                                                       3586, 9931, 973, 1001,
                                                                       4384, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 11569, 0, 3,
                                                                       10057, 3712, 10225, 1057,
                                                                       1093, 4468, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 11785, 0, 3,
                                                                       10225, 3796, 10393, 1093,
                                                                       1129, 4576, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12001, 0, 3,
                                                                       10393, 3880, 10561, 1129,
                                                                       1165, 4684, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12217, 0, 3,
                                                                       10561, 3964, 10729, 1165,
                                                                       1201, 4792, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12433, 0, 3,
                                                                       10729, 4048, 10897, 1201,
                                                                       1237, 4900, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12649, 0, 3,
                                                                       10897, 4132, 11065, 1237,
                                                                       1273, 5008, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12865, 0, 3,
                                                                       11065, 4216, 11233, 1273,
                                                                       1309, 5116, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 13081, 0, 3,
                                                                       11233, 4300, 11401, 1309,
                                                                       1345, 5224, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 13297, 0, 3,
                                                                       11569, 4468, 11785, 1417,
                                                                       1462, 5332, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 13567, 0, 3,
                                                                       11785, 4576, 12001, 1462,
                                                                       1507, 5467, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 13837, 0, 3,
                                                                       12001, 4684, 12217, 1507,
                                                                       1552, 5602, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 14107, 0, 3,
                                                                       12217, 4792, 12433, 1552,
                                                                       1597, 5737, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 14377, 0, 3,
                                                                       12433, 4900, 12649, 1597,
                                                                       1642, 5872, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 14647, 0, 3,
                                                                       12649, 5008, 12865, 1642,
                                                                       1687, 6007, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 14917, 0, 3,
                                                                       12865, 5116, 13081, 1687,
                                                                       1732, 6142, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15187, 3, 1822,
                                                                       1825, 6289, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15197, 3, 1825,
                                                                       1828, 6295, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15207, 3, 1828,
                                                                       1831, 6301, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15217, 3, 1831,
                                                                       1834, 6307, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15227, 3, 1834,
                                                                       1837, 6313, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15237, 3, 1837,
                                                                       1840, 6319, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15247, 3, 1840,
                                                                       1843, 6325, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15257, 3, 1843,
                                                                       1846, 6331, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15267, 3, 1846,
                                                                       1849, 6337, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15277, 3, 1849,
                                                                       1852, 6343, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15287, 3, 1852,
                                                                       1855, 6349, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15297, 3, 1855,
                                                                       1858, 6355, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15307, 3, 1858,
                                                                       1861, 6361, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15317, 0, 3,
                                                                       15187, 6289, 15197, 1867,
                                                                       1876, 6403, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15347, 0, 3,
                                                                       15197, 6295, 15207, 1876,
                                                                       1885, 6421, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15377, 0, 3,
                                                                       15207, 6301, 15217, 1885,
                                                                       1894, 6439, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15407, 0, 3,
                                                                       15217, 6307, 15227, 1894,
                                                                       1903, 6457, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15437, 0, 3,
                                                                       15227, 6313, 15237, 1903,
                                                                       1912, 6475, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15467, 0, 3,
                                                                       15237, 6319, 15247, 1912,
                                                                       1921, 6493, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15497, 0, 3,
                                                                       15247, 6325, 15257, 1921,
                                                                       1930, 6511, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15527, 0, 3,
                                                                       15257, 6331, 15267, 1930,
                                                                       1939, 6529, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15557, 0, 3,
                                                                       15267, 6337, 15277, 1939,
                                                                       1948, 6547, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15587, 0, 3,
                                                                       15277, 6343, 15287, 1948,
                                                                       1957, 6565, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15617, 0, 3,
                                                                       15287, 6349, 15297, 1957,
                                                                       1966, 6583, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15647, 0, 3,
                                                                       15297, 6355, 15307, 1966,
                                                                       1975, 6601, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 15677, 0, 3,
                                                                       15317, 6403, 15347, 1993,
                                                                       2011, 6691, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 15737, 0, 3,
                                                                       15347, 6421, 15377, 2011,
                                                                       2029, 6727, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 15797, 0, 3,
                                                                       15377, 6439, 15407, 2029,
                                                                       2047, 6763, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 15857, 0, 3,
                                                                       15407, 6457, 15437, 2047,
                                                                       2065, 6799, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 15917, 0, 3,
                                                                       15437, 6475, 15467, 2065,
                                                                       2083, 6835, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 15977, 0, 3,
                                                                       15467, 6493, 15497, 2083,
                                                                       2101, 6871, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16037, 0, 3,
                                                                       15497, 6511, 15527, 2101,
                                                                       2119, 6907, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16097, 0, 3,
                                                                       15527, 6529, 15557, 2119,
                                                                       2137, 6943, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16157, 0, 3,
                                                                       15557, 6547, 15587, 2137,
                                                                       2155, 6979, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16217, 0, 3,
                                                                       15587, 6565, 15617, 2155,
                                                                       2173, 7015, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16277, 0, 3,
                                                                       15617, 6583, 15647, 2173,
                                                                       2191, 7051, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 16337, 0, 3,
                                                                       15677, 6691, 15737, 2227,
                                                                       2257, 7207, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 16437, 0, 3,
                                                                       15737, 6727, 15797, 2257,
                                                                       2287, 7267, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 16537, 0, 3,
                                                                       15797, 6763, 15857, 2287,
                                                                       2317, 7327, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 16637, 0, 3,
                                                                       15857, 6799, 15917, 2317,
                                                                       2347, 7387, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 16737, 0, 3,
                                                                       15917, 6835, 15977, 2347,
                                                                       2377, 7447, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 16837, 0, 3,
                                                                       15977, 6871, 16037, 2377,
                                                                       2407, 7507, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 16937, 0, 3,
                                                                       16037, 6907, 16097, 2407,
                                                                       2437, 7567, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17037, 0, 3,
                                                                       16097, 6943, 16157, 2437,
                                                                       2467, 7627, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17137, 0, 3,
                                                                       16157, 6979, 16217, 2467,
                                                                       2497, 7687, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17237, 0, 3,
                                                                       16217, 7015, 16277, 2497,
                                                                       2527, 7747, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17337, 0, 3,
                                                                       16337, 7207, 16437, 2587,
                                                                       2632, 7987, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17487, 0, 3,
                                                                       16437, 7267, 16537, 2632,
                                                                       2677, 8077, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17637, 0, 3,
                                                                       16537, 7327, 16637, 2677,
                                                                       2722, 8167, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17787, 0, 3,
                                                                       16637, 7387, 16737, 2722,
                                                                       2767, 8257, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17937, 0, 3,
                                                                       16737, 7447, 16837, 2767,
                                                                       2812, 8347, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18087, 0, 3,
                                                                       16837, 7507, 16937, 2812,
                                                                       2857, 8437, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18237, 0, 3,
                                                                       16937, 7567, 17037, 2857,
                                                                       2902, 8527, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18387, 0, 3,
                                                                       17037, 7627, 17137, 2902,
                                                                       2947, 8617, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18537, 0, 3,
                                                                       17137, 7687, 17237, 2947,
                                                                       2992, 8707, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 18687, 0, 3,
                                                                       17337, 7987, 17487, 3082,
                                                                       3145, 9049, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 18897, 0, 3,
                                                                       17487, 8077, 17637, 3145,
                                                                       3208, 9175, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19107, 0, 3,
                                                                       17637, 8167, 17787, 3208,
                                                                       3271, 9301, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19317, 0, 3,
                                                                       17787, 8257, 17937, 3271,
                                                                       3334, 9427, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19527, 0, 3,
                                                                       17937, 8347, 18087, 3334,
                                                                       3397, 9553, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19737, 0, 3,
                                                                       18087, 8437, 18237, 3397,
                                                                       3460, 9679, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19947, 0, 3,
                                                                       18237, 8527, 18387, 3460,
                                                                       3523, 9805, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 20157, 0, 3,
                                                                       18387, 8617, 18537, 3523,
                                                                       3586, 9931, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 20367, 0, 3,
                                                                       18687, 9049, 18897, 3712,
                                                                       3796, 10393, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 20647, 0, 3,
                                                                       18897, 9175, 19107, 3796,
                                                                       3880, 10561, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 20927, 0, 3,
                                                                       19107, 9301, 19317, 3880,
                                                                       3964, 10729, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 21207, 0, 3,
                                                                       19317, 9427, 19527, 3964,
                                                                       4048, 10897, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 21487, 0, 3,
                                                                       19527, 9553, 19737, 4048,
                                                                       4132, 11065, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 21767, 0, 3,
                                                                       19737, 9679, 19947, 4132,
                                                                       4216, 11233, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 22047, 0, 3,
                                                                       19947, 9805, 20157, 4216,
                                                                       4300, 11401, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 22327, 0, 3,
                                                                       20367, 10393, 20647, 4468,
                                                                       4576, 12001, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 22687, 0, 3,
                                                                       20647, 10561, 20927, 4576,
                                                                       4684, 12217, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 23047, 0, 3,
                                                                       20927, 10729, 21207, 4684,
                                                                       4792, 12433, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 23407, 0, 3,
                                                                       21207, 10897, 21487, 4792,
                                                                       4900, 12649, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 23767, 0, 3,
                                                                       21487, 11065, 21767, 4900,
                                                                       5008, 12865, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 24127, 0, 3,
                                                                       21767, 11233, 22047, 5008,
                                                                       5116, 13081, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 24487, 0, 3,
                                                                       22327, 12001, 22687, 5332,
                                                                       5467, 13837, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 24937, 0, 3,
                                                                       22687, 12217, 23047, 5467,
                                                                       5602, 14107, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 25387, 0, 3,
                                                                       23047, 12433, 23407, 5602,
                                                                       5737, 14377, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 25837, 0, 3,
                                                                       23407, 12649, 23767, 5737,
                                                                       5872, 14647, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 26287, 0, 3,
                                                                       23767, 12865, 24127, 5872,
                                                                       6007, 14917, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26737, 3, 6277,
                                                                       6283, 15187, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26752, 3, 6283,
                                                                       6289, 15197, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26767, 3, 6289,
                                                                       6295, 15207, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26782, 3, 6295,
                                                                       6301, 15217, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26797, 3, 6301,
                                                                       6307, 15227, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26812, 3, 6307,
                                                                       6313, 15237, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26827, 3, 6313,
                                                                       6319, 15247, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26842, 3, 6319,
                                                                       6325, 15257, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26857, 3, 6325,
                                                                       6331, 15267, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26872, 3, 6331,
                                                                       6337, 15277, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26887, 3, 6337,
                                                                       6343, 15287, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26902, 3, 6343,
                                                                       6349, 15297, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26917, 3, 6349,
                                                                       6355, 15307, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26932, 0, 3,
                                                                       26737, 15187, 26752, 6367,
                                                                       6385, 15317, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26977, 0, 3,
                                                                       26752, 15197, 26767, 6385,
                                                                       6403, 15347, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27022, 0, 3,
                                                                       26767, 15207, 26782, 6403,
                                                                       6421, 15377, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27067, 0, 3,
                                                                       26782, 15217, 26797, 6421,
                                                                       6439, 15407, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27112, 0, 3,
                                                                       26797, 15227, 26812, 6439,
                                                                       6457, 15437, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27157, 0, 3,
                                                                       26812, 15237, 26827, 6457,
                                                                       6475, 15467, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27202, 0, 3,
                                                                       26827, 15247, 26842, 6475,
                                                                       6493, 15497, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27247, 0, 3,
                                                                       26842, 15257, 26857, 6493,
                                                                       6511, 15527, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27292, 0, 3,
                                                                       26857, 15267, 26872, 6511,
                                                                       6529, 15557, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27337, 0, 3,
                                                                       26872, 15277, 26887, 6529,
                                                                       6547, 15587, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27382, 0, 3,
                                                                       26887, 15287, 26902, 6547,
                                                                       6565, 15617, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27427, 0, 3,
                                                                       26902, 15297, 26917, 6565,
                                                                       6583, 15647, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 27472, 0, 3,
                                                                       26932, 15317, 26977, 6619,
                                                                       6655, 15677, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 27562, 0, 3,
                                                                       26977, 15347, 27022, 6655,
                                                                       6691, 15737, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 27652, 0, 3,
                                                                       27022, 15377, 27067, 6691,
                                                                       6727, 15797, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 27742, 0, 3,
                                                                       27067, 15407, 27112, 6727,
                                                                       6763, 15857, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 27832, 0, 3,
                                                                       27112, 15437, 27157, 6763,
                                                                       6799, 15917, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 27922, 0, 3,
                                                                       27157, 15467, 27202, 6799,
                                                                       6835, 15977, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 28012, 0, 3,
                                                                       27202, 15497, 27247, 6835,
                                                                       6871, 16037, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 28102, 0, 3,
                                                                       27247, 15527, 27292, 6871,
                                                                       6907, 16097, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 28192, 0, 3,
                                                                       27292, 15557, 27337, 6907,
                                                                       6943, 16157, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 28282, 0, 3,
                                                                       27337, 15587, 27382, 6943,
                                                                       6979, 16217, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 28372, 0, 3,
                                                                       27382, 15617, 27427, 6979,
                                                                       7015, 16277, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 28462, 0, 3,
                                                                       27472, 15677, 27562, 7087,
                                                                       7147, 16337, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 28612, 0, 3,
                                                                       27562, 15737, 27652, 7147,
                                                                       7207, 16437, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 28762, 0, 3,
                                                                       27652, 15797, 27742, 7207,
                                                                       7267, 16537, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 28912, 0, 3,
                                                                       27742, 15857, 27832, 7267,
                                                                       7327, 16637, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 29062, 0, 3,
                                                                       27832, 15917, 27922, 7327,
                                                                       7387, 16737, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 29212, 0, 3,
                                                                       27922, 15977, 28012, 7387,
                                                                       7447, 16837, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 29362, 0, 3,
                                                                       28012, 16037, 28102, 7447,
                                                                       7507, 16937, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 29512, 0, 3,
                                                                       28102, 16097, 28192, 7507,
                                                                       7567, 17037, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 29662, 0, 3,
                                                                       28192, 16157, 28282, 7567,
                                                                       7627, 17137, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 29812, 0, 3,
                                                                       28282, 16217, 28372, 7627,
                                                                       7687, 17237, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 29962, 0, 3,
                                                                       28462, 16337, 28612, 7807,
                                                                       7897, 17337, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 30187, 0, 3,
                                                                       28612, 16437, 28762, 7897,
                                                                       7987, 17487, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 30412, 0, 3,
                                                                       28762, 16537, 28912, 7987,
                                                                       8077, 17637, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 30637, 0, 3,
                                                                       28912, 16637, 29062, 8077,
                                                                       8167, 17787, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 30862, 0, 3,
                                                                       29062, 16737, 29212, 8167,
                                                                       8257, 17937, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 31087, 0, 3,
                                                                       29212, 16837, 29362, 8257,
                                                                       8347, 18087, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 31312, 0, 3,
                                                                       29362, 16937, 29512, 8347,
                                                                       8437, 18237, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 31537, 0, 3,
                                                                       29512, 17037, 29662, 8437,
                                                                       8527, 18387, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 31762, 0, 3,
                                                                       29662, 17137, 29812, 8527,
                                                                       8617, 18537, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 31987, 0, 3,
                                                                       29962, 17337, 30187, 8797,
                                                                       8923, 18687, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 32302, 0, 3,
                                                                       30187, 17487, 30412, 8923,
                                                                       9049, 18897, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 32617, 0, 3,
                                                                       30412, 17637, 30637, 9049,
                                                                       9175, 19107, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 32932, 0, 3,
                                                                       30637, 17787, 30862, 9175,
                                                                       9301, 19317, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 33247, 0, 3,
                                                                       30862, 17937, 31087, 9301,
                                                                       9427, 19527, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 33562, 0, 3,
                                                                       31087, 18087, 31312, 9427,
                                                                       9553, 19737, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 33877, 0, 3,
                                                                       31312, 18237, 31537, 9553,
                                                                       9679, 19947, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 34192, 0, 3,
                                                                       31537, 18387, 31762, 9679,
                                                                       9805, 20157, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 34507, 0, 3,
                                                                       31987, 18687, 32302,
                                                                       10057, 10225, 20367,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 34927, 0, 3,
                                                                       32302, 18897, 32617,
                                                                       10225, 10393, 20647,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 35347, 0, 3,
                                                                       32617, 19107, 32932,
                                                                       10393, 10561, 20927,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 35767, 0, 3,
                                                                       32932, 19317, 33247,
                                                                       10561, 10729, 21207,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 36187, 0, 3,
                                                                       33247, 19527, 33562,
                                                                       10729, 10897, 21487,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 36607, 0, 3,
                                                                       33562, 19737, 33877,
                                                                       10897, 11065, 21767,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 37027, 0, 3,
                                                                       33877, 19947, 34192,
                                                                       11065, 11233, 22047,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 37447, 0, 3,
                                                                       34507, 20367, 34927,
                                                                       11569, 11785, 22327,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 37987, 0, 3,
                                                                       34927, 20647, 35347,
                                                                       11785, 12001, 22687,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 38527, 0, 3,
                                                                       35347, 20927, 35767,
                                                                       12001, 12217, 23047,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 39067, 0, 3,
                                                                       35767, 21207, 36187,
                                                                       12217, 12433, 23407,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 39607, 0, 3,
                                                                       36187, 21487, 36607,
                                                                       12433, 12649, 23767,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 40147, 0, 3,
                                                                       36607, 21767, 37027,
                                                                       12649, 12865, 24127,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 40687, 0, 3,
                                                                       37447, 22327, 37987,
                                                                       13297, 13567, 24487,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 41362, 0, 3,
                                                                       37987, 22687, 38527,
                                                                       13567, 13837, 24937,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 42037, 0, 3,
                                                                       38527, 23047, 39067,
                                                                       13837, 14107, 25387,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 42712, 0, 3,
                                                                       39067, 23407, 39607,
                                                                       14107, 14377, 25837,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 43387, 0, 3,
                                                                       39607, 23767, 40147,
                                                                       14377, 14647, 26287,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44062, 3, 15187,
                                                                       15197, 26767, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44083, 3, 15197,
                                                                       15207, 26782, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44104, 3, 15207,
                                                                       15217, 26797, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44125, 3, 15217,
                                                                       15227, 26812, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44146, 3, 15227,
                                                                       15237, 26827, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44167, 3, 15237,
                                                                       15247, 26842, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44188, 3, 15247,
                                                                       15257, 26857, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44209, 3, 15257,
                                                                       15267, 26872, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44230, 3, 15267,
                                                                       15277, 26887, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44251, 3, 15277,
                                                                       15287, 26902, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44272, 3, 15287,
                                                                       15297, 26917, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 44293, 0, 3,
                                                                       44062, 26767, 44083,
                                                                       15317, 15347, 27022,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 44356, 0, 3,
                                                                       44083, 26782, 44104,
                                                                       15347, 15377, 27067,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 44419, 0, 3,
                                                                       44104, 26797, 44125,
                                                                       15377, 15407, 27112,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 44482, 0, 3,
                                                                       44125, 26812, 44146,
                                                                       15407, 15437, 27157,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 44545, 0, 3,
                                                                       44146, 26827, 44167,
                                                                       15437, 15467, 27202,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 44608, 0, 3,
                                                                       44167, 26842, 44188,
                                                                       15467, 15497, 27247,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 44671, 0, 3,
                                                                       44188, 26857, 44209,
                                                                       15497, 15527, 27292,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 44734, 0, 3,
                                                                       44209, 26872, 44230,
                                                                       15527, 15557, 27337,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 44797, 0, 3,
                                                                       44230, 26887, 44251,
                                                                       15557, 15587, 27382,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 44860, 0, 3,
                                                                       44251, 26902, 44272,
                                                                       15587, 15617, 27427,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 44923, 0, 3,
                                                                       44293, 27022, 44356,
                                                                       15677, 15737, 27652,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 45049, 0, 3,
                                                                       44356, 27067, 44419,
                                                                       15737, 15797, 27742,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 45175, 0, 3,
                                                                       44419, 27112, 44482,
                                                                       15797, 15857, 27832,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 45301, 0, 3,
                                                                       44482, 27157, 44545,
                                                                       15857, 15917, 27922,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 45427, 0, 3,
                                                                       44545, 27202, 44608,
                                                                       15917, 15977, 28012,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 45553, 0, 3,
                                                                       44608, 27247, 44671,
                                                                       15977, 16037, 28102,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 45679, 0, 3,
                                                                       44671, 27292, 44734,
                                                                       16037, 16097, 28192,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 45805, 0, 3,
                                                                       44734, 27337, 44797,
                                                                       16097, 16157, 28282,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 45931, 0, 3,
                                                                       44797, 27382, 44860,
                                                                       16157, 16217, 28372,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 46057, 0, 3,
                                                                       44923, 27652, 45049,
                                                                       16337, 16437, 28762,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 46267, 0, 3,
                                                                       45049, 27742, 45175,
                                                                       16437, 16537, 28912,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 46477, 0, 3,
                                                                       45175, 27832, 45301,
                                                                       16537, 16637, 29062,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 46687, 0, 3,
                                                                       45301, 27922, 45427,
                                                                       16637, 16737, 29212,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 46897, 0, 3,
                                                                       45427, 28012, 45553,
                                                                       16737, 16837, 29362,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 47107, 0, 3,
                                                                       45553, 28102, 45679,
                                                                       16837, 16937, 29512,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 47317, 0, 3,
                                                                       45679, 28192, 45805,
                                                                       16937, 17037, 29662,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 47527, 0, 3,
                                                                       45805, 28282, 45931,
                                                                       17037, 17137, 29812,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 47737, 0, 3,
                                                                       46057, 28762, 46267,
                                                                       17337, 17487, 30412,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 48052, 0, 3,
                                                                       46267, 28912, 46477,
                                                                       17487, 17637, 30637,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 48367, 0, 3,
                                                                       46477, 29062, 46687,
                                                                       17637, 17787, 30862,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 48682, 0, 3,
                                                                       46687, 29212, 46897,
                                                                       17787, 17937, 31087,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 48997, 0, 3,
                                                                       46897, 29362, 47107,
                                                                       17937, 18087, 31312,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 49312, 0, 3,
                                                                       47107, 29512, 47317,
                                                                       18087, 18237, 31537,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 49627, 0, 3,
                                                                       47317, 29662, 47527,
                                                                       18237, 18387, 31762,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 49942, 0, 3,
                                                                       47737, 30412, 48052,
                                                                       18687, 18897, 32617,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 50383, 0, 3,
                                                                       48052, 30637, 48367,
                                                                       18897, 19107, 32932,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 50824, 0, 3,
                                                                       48367, 30862, 48682,
                                                                       19107, 19317, 33247,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 51265, 0, 3,
                                                                       48682, 31087, 48997,
                                                                       19317, 19527, 33562,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 51706, 0, 3,
                                                                       48997, 31312, 49312,
                                                                       19527, 19737, 33877,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 52147, 0, 3,
                                                                       49312, 31537, 49627,
                                                                       19737, 19947, 34192,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 52588, 0, 3,
                                                                       49942, 32617, 50383,
                                                                       20367, 20647, 35347,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 53176, 0, 3,
                                                                       50383, 32932, 50824,
                                                                       20647, 20927, 35767,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 53764, 0, 3,
                                                                       50824, 33247, 51265,
                                                                       20927, 21207, 36187,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 54352, 0, 3,
                                                                       51265, 33562, 51706,
                                                                       21207, 21487, 36607,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 54940, 0, 3,
                                                                       51706, 33877, 52147,
                                                                       21487, 21767, 37027,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 55528, 0, 3,
                                                                       52588, 35347, 53176,
                                                                       22327, 22687, 38527,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 56284, 0, 3,
                                                                       53176, 35767, 53764,
                                                                       22687, 23047, 39067,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 57040, 0, 3,
                                                                       53764, 36187, 54352,
                                                                       23047, 23407, 39607,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 57796, 0, 3,
                                                                       54352, 36607, 54940,
                                                                       23407, 23767, 40147,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 58552, 0, 3,
                                                                       55528, 38527, 56284,
                                                                       24487, 24937, 42037,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 59497, 0, 3,
                                                                       56284, 39067, 57040,
                                                                       24937, 25387, 42712,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 60442, 0, 3,
                                                                       57040, 39607, 57796,
                                                                       25387, 25837, 43387,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61387, 3, 26737,
                                                                       26752, 44062, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61415, 3, 26752,
                                                                       26767, 44083, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61443, 3, 26767,
                                                                       26782, 44104, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61471, 3, 26782,
                                                                       26797, 44125, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61499, 3, 26797,
                                                                       26812, 44146, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61527, 3, 26812,
                                                                       26827, 44167, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61555, 3, 26827,
                                                                       26842, 44188, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61583, 3, 26842,
                                                                       26857, 44209, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61611, 3, 26857,
                                                                       26872, 44230, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61639, 3, 26872,
                                                                       26887, 44251, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61667, 3, 26887,
                                                                       26902, 44272, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 61695, 0, 3,
                                                                       61387, 44062, 61415,
                                                                       26932, 26977, 44293,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 61779, 0, 3,
                                                                       61415, 44083, 61443,
                                                                       26977, 27022, 44356,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 61863, 0, 3,
                                                                       61443, 44104, 61471,
                                                                       27022, 27067, 44419,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 61947, 0, 3,
                                                                       61471, 44125, 61499,
                                                                       27067, 27112, 44482,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62031, 0, 3,
                                                                       61499, 44146, 61527,
                                                                       27112, 27157, 44545,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62115, 0, 3,
                                                                       61527, 44167, 61555,
                                                                       27157, 27202, 44608,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62199, 0, 3,
                                                                       61555, 44188, 61583,
                                                                       27202, 27247, 44671,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62283, 0, 3,
                                                                       61583, 44209, 61611,
                                                                       27247, 27292, 44734,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62367, 0, 3,
                                                                       61611, 44230, 61639,
                                                                       27292, 27337, 44797,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62451, 0, 3,
                                                                       61639, 44251, 61667,
                                                                       27337, 27382, 44860,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 62535, 0, 3,
                                                                       61695, 44293, 61779,
                                                                       27472, 27562, 44923,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 62703, 0, 3,
                                                                       61779, 44356, 61863,
                                                                       27562, 27652, 45049,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 62871, 0, 3,
                                                                       61863, 44419, 61947,
                                                                       27652, 27742, 45175,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 63039, 0, 3,
                                                                       61947, 44482, 62031,
                                                                       27742, 27832, 45301,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 63207, 0, 3,
                                                                       62031, 44545, 62115,
                                                                       27832, 27922, 45427,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 63375, 0, 3,
                                                                       62115, 44608, 62199,
                                                                       27922, 28012, 45553,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 63543, 0, 3,
                                                                       62199, 44671, 62283,
                                                                       28012, 28102, 45679,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 63711, 0, 3,
                                                                       62283, 44734, 62367,
                                                                       28102, 28192, 45805,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 63879, 0, 3,
                                                                       62367, 44797, 62451,
                                                                       28192, 28282, 45931,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 64047, 0, 3,
                                                                       62535, 44923, 62703,
                                                                       28462, 28612, 46057,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 64327, 0, 3,
                                                                       62703, 45049, 62871,
                                                                       28612, 28762, 46267,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 64607, 0, 3,
                                                                       62871, 45175, 63039,
                                                                       28762, 28912, 46477,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 64887, 0, 3,
                                                                       63039, 45301, 63207,
                                                                       28912, 29062, 46687,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 65167, 0, 3,
                                                                       63207, 45427, 63375,
                                                                       29062, 29212, 46897,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 65447, 0, 3,
                                                                       63375, 45553, 63543,
                                                                       29212, 29362, 47107,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 65727, 0, 3,
                                                                       63543, 45679, 63711,
                                                                       29362, 29512, 47317,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 66007, 0, 3,
                                                                       63711, 45805, 63879,
                                                                       29512, 29662, 47527,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 66287, 0, 3,
                                                                       64047, 46057, 64327,
                                                                       29962, 30187, 47737,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 66707, 0, 3,
                                                                       64327, 46267, 64607,
                                                                       30187, 30412, 48052,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 67127, 0, 3,
                                                                       64607, 46477, 64887,
                                                                       30412, 30637, 48367,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 67547, 0, 3,
                                                                       64887, 46687, 65167,
                                                                       30637, 30862, 48682,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 67967, 0, 3,
                                                                       65167, 46897, 65447,
                                                                       30862, 31087, 48997,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 68387, 0, 3,
                                                                       65447, 47107, 65727,
                                                                       31087, 31312, 49312,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 68807, 0, 3,
                                                                       65727, 47317, 66007,
                                                                       31312, 31537, 49627,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 69227, 0, 3,
                                                                       66287, 47737, 66707,
                                                                       31987, 32302, 49942,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 69815, 0, 3,
                                                                       66707, 48052, 67127,
                                                                       32302, 32617, 50383,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 70403, 0, 3,
                                                                       67127, 48367, 67547,
                                                                       32617, 32932, 50824,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 70991, 0, 3,
                                                                       67547, 48682, 67967,
                                                                       32932, 33247, 51265,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 71579, 0, 3,
                                                                       67967, 48997, 68387,
                                                                       33247, 33562, 51706,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 72167, 0, 3,
                                                                       68387, 49312, 68807,
                                                                       33562, 33877, 52147,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 72755, 0, 3,
                                                                       69227, 49942, 69815,
                                                                       34507, 34927, 52588,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 73539, 0, 3,
                                                                       69815, 50383, 70403,
                                                                       34927, 35347, 53176,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 74323, 0, 3,
                                                                       70403, 50824, 70991,
                                                                       35347, 35767, 53764,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 75107, 0, 3,
                                                                       70991, 51265, 71579,
                                                                       35767, 36187, 54352,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 75891, 0, 3,
                                                                       71579, 51706, 72167,
                                                                       36187, 36607, 54940,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 76675, 0, 3,
                                                                       72755, 52588, 73539,
                                                                       37447, 37987, 55528,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 77683, 0, 3,
                                                                       73539, 53176, 74323,
                                                                       37987, 38527, 56284,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 78691, 0, 3,
                                                                       74323, 53764, 75107,
                                                                       38527, 39067, 57040,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 79699, 0, 3,
                                                                       75107, 54352, 75891,
                                                                       39067, 39607, 57796,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 80707, 0, 3,
                                                                       76675, 55528, 77683,
                                                                       40687, 41362, 58552,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 81967, 0, 3,
                                                                       77683, 56284, 78691,
                                                                       41362, 42037, 59497,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 83227, 0, 3,
                                                                       78691, 57040, 79699,
                                                                       42037, 42712, 60442,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 84487, 3, 44062,
                                                                       44083, 61443, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 84523, 3, 44083,
                                                                       44104, 61471, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 84559, 3, 44104,
                                                                       44125, 61499, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 84595, 3, 44125,
                                                                       44146, 61527, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 84631, 3, 44146,
                                                                       44167, 61555, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 84667, 3, 44167,
                                                                       44188, 61583, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 84703, 3, 44188,
                                                                       44209, 61611, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 84739, 3, 44209,
                                                                       44230, 61639, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 84775, 3, 44230,
                                                                       44251, 61667, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 84811, 0, 3,
                                                                       84487, 61443, 84523,
                                                                       44293, 44356, 61863,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 84919, 0, 3,
                                                                       84523, 61471, 84559,
                                                                       44356, 44419, 61947,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 85027, 0, 3,
                                                                       84559, 61499, 84595,
                                                                       44419, 44482, 62031,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 85135, 0, 3,
                                                                       84595, 61527, 84631,
                                                                       44482, 44545, 62115,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 85243, 0, 3,
                                                                       84631, 61555, 84667,
                                                                       44545, 44608, 62199,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 85351, 0, 3,
                                                                       84667, 61583, 84703,
                                                                       44608, 44671, 62283,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 85459, 0, 3,
                                                                       84703, 61611, 84739,
                                                                       44671, 44734, 62367,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 85567, 0, 3,
                                                                       84739, 61639, 84775,
                                                                       44734, 44797, 62451,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 85675, 0, 3,
                                                                       84811, 61863, 84919,
                                                                       44923, 45049, 62871,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 85891, 0, 3,
                                                                       84919, 61947, 85027,
                                                                       45049, 45175, 63039,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 86107, 0, 3,
                                                                       85027, 62031, 85135,
                                                                       45175, 45301, 63207,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 86323, 0, 3,
                                                                       85135, 62115, 85243,
                                                                       45301, 45427, 63375,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 86539, 0, 3,
                                                                       85243, 62199, 85351,
                                                                       45427, 45553, 63543,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 86755, 0, 3,
                                                                       85351, 62283, 85459,
                                                                       45553, 45679, 63711,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 86971, 0, 3,
                                                                       85459, 62367, 85567,
                                                                       45679, 45805, 63879,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 87187, 0, 3,
                                                                       85675, 62871, 85891,
                                                                       46057, 46267, 64607,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 87547, 0, 3,
                                                                       85891, 63039, 86107,
                                                                       46267, 46477, 64887,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 87907, 0, 3,
                                                                       86107, 63207, 86323,
                                                                       46477, 46687, 65167,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 88267, 0, 3,
                                                                       86323, 63375, 86539,
                                                                       46687, 46897, 65447,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 88627, 0, 3,
                                                                       86539, 63543, 86755,
                                                                       46897, 47107, 65727,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 88987, 0, 3,
                                                                       86755, 63711, 86971,
                                                                       47107, 47317, 66007,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 89347, 0, 3,
                                                                       87187, 64607, 87547,
                                                                       47737, 48052, 67127,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 89887, 0, 3,
                                                                       87547, 64887, 87907,
                                                                       48052, 48367, 67547,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 90427, 0, 3,
                                                                       87907, 65167, 88267,
                                                                       48367, 48682, 67967,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 90967, 0, 3,
                                                                       88267, 65447, 88627,
                                                                       48682, 48997, 68387,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 91507, 0, 3,
                                                                       88627, 65727, 88987,
                                                                       48997, 49312, 68807,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 92047, 0, 3,
                                                                       89347, 67127, 89887,
                                                                       49942, 50383, 70403,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 92803, 0, 3,
                                                                       89887, 67547, 90427,
                                                                       50383, 50824, 70991,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 93559, 0, 3,
                                                                       90427, 67967, 90967,
                                                                       50824, 51265, 71579,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 94315, 0, 3,
                                                                       90967, 68387, 91507,
                                                                       51265, 51706, 72167,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 95071, 0, 3,
                                                                       92047, 70403, 92803,
                                                                       52588, 53176, 74323,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 96079, 0, 3,
                                                                       92803, 70991, 93559,
                                                                       53176, 53764, 75107,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 97087, 0, 3,
                                                                       93559, 71579, 94315,
                                                                       53764, 54352, 75891,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 98095, 0, 3,
                                                                       95071, 74323, 96079,
                                                                       55528, 56284, 78691,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 99391, 0, 3,
                                                                       96079, 75107, 97087,
                                                                       56284, 57040, 79699,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 100687, 0, 3,
                                                                       98095, 78691, 99391,
                                                                       58552, 59497, 83227,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 102307, 3, 61387,
                                                                       61415, 84487, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 102352, 3, 61415,
                                                                       61443, 84523, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 102397, 3, 61443,
                                                                       61471, 84559, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 102442, 3, 61471,
                                                                       61499, 84595, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 102487, 3, 61499,
                                                                       61527, 84631, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 102532, 3, 61527,
                                                                       61555, 84667, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 102577, 3, 61555,
                                                                       61583, 84703, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 102622, 3, 61583,
                                                                       61611, 84739, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 102667, 3, 61611,
                                                                       61639, 84775, ncols,
                                                                       gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 102712, 0, 3,
                                                                       102307, 84487, 102352,
                                                                       61695, 61779, 84811,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 102847, 0, 3,
                                                                       102352, 84523, 102397,
                                                                       61779, 61863, 84919,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 102982, 0, 3,
                                                                       102397, 84559, 102442,
                                                                       61863, 61947, 85027,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 103117, 0, 3,
                                                                       102442, 84595, 102487,
                                                                       61947, 62031, 85135,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 103252, 0, 3,
                                                                       102487, 84631, 102532,
                                                                       62031, 62115, 85243,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 103387, 0, 3,
                                                                       102532, 84667, 102577,
                                                                       62115, 62199, 85351,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 103522, 0, 3,
                                                                       102577, 84703, 102622,
                                                                       62199, 62283, 85459,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 103657, 0, 3,
                                                                       102622, 84739, 102667,
                                                                       62283, 62367, 85567,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 103792, 0, 3,
                                                                       102712, 84811, 102847,
                                                                       62535, 62703, 85675,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 104062, 0, 3,
                                                                       102847, 84919, 102982,
                                                                       62703, 62871, 85891,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 104332, 0, 3,
                                                                       102982, 85027, 103117,
                                                                       62871, 63039, 86107,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 104602, 0, 3,
                                                                       103117, 85135, 103252,
                                                                       63039, 63207, 86323,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 104872, 0, 3,
                                                                       103252, 85243, 103387,
                                                                       63207, 63375, 86539,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 105142, 0, 3,
                                                                       103387, 85351, 103522,
                                                                       63375, 63543, 86755,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 105412, 0, 3,
                                                                       103522, 85459, 103657,
                                                                       63543, 63711, 86971,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 105682, 0, 3,
                                                                       103792, 85675, 104062,
                                                                       64047, 64327, 87187,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 106132, 0, 3,
                                                                       104062, 85891, 104332,
                                                                       64327, 64607, 87547,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 106582, 0, 3,
                                                                       104332, 86107, 104602,
                                                                       64607, 64887, 87907,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 107032, 0, 3,
                                                                       104602, 86323, 104872,
                                                                       64887, 65167, 88267,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 107482, 0, 3,
                                                                       104872, 86539, 105142,
                                                                       65167, 65447, 88627,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 107932, 0, 3,
                                                                       105142, 86755, 105412,
                                                                       65447, 65727, 88987,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 108382, 0, 3,
                                                                       105682, 87187, 106132,
                                                                       66287, 66707, 89347,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 109057, 0, 3,
                                                                       106132, 87547, 106582,
                                                                       66707, 67127, 89887,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 109732, 0, 3,
                                                                       106582, 87907, 107032,
                                                                       67127, 67547, 90427,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 110407, 0, 3,
                                                                       107032, 88267, 107482,
                                                                       67547, 67967, 90967,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 111082, 0, 3,
                                                                       107482, 88627, 107932,
                                                                       67967, 68387, 91507,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 111757, 0, 3,
                                                                       108382, 89347, 109057,
                                                                       69227, 69815, 92047,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 112702, 0, 3,
                                                                       109057, 89887, 109732,
                                                                       69815, 70403, 92803,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 113647, 0, 3,
                                                                       109732, 90427, 110407,
                                                                       70403, 70991, 93559,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 114592, 0, 3,
                                                                       110407, 90967, 111082,
                                                                       70991, 71579, 94315,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 115537, 0, 3,
                                                                       111757, 92047, 112702,
                                                                       72755, 73539, 95071,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 116797, 0, 3,
                                                                       112702, 92803, 113647,
                                                                       73539, 74323, 96079,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 118057, 0, 3,
                                                                       113647, 93559, 114592,
                                                                       74323, 75107, 97087,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 119317, 0, 3,
                                                                       115537, 95071, 116797,
                                                                       76675, 77683, 98095,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 120937, 0, 3,
                                                                       116797, 96079, 118057,
                                                                       77683, 78691, 99391,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 122557, 0, 3,
                                                                       119317, 98095, 120937,
                                                                       80707, 81967, 100687,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 124582, 115537, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 126318, 119317, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 128550, 122557, 2025, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 125842, 124582, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 127938, 126318, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 130575, 128550, 45, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 131340, 125842, 127938, 17,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 132768, 127938, 130575, 17,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 134604, 131340, 132768, 17,
                                             nmax);

        simdtrf::transform_d_inner(buffer, 137460, 134604, 28, 17, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 137460, 85, nmax);
    }

    for (size_t m = 0; m < 1105; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
