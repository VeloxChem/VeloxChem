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


#include "SimdThreeCenterElectronRepulsionRsRecIFL.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
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
#include "SimdTransferIF.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_ifl_three_center_electron_repulsion(double               *values,
                                               const size_t          npairs,
                                               const size_t          natoms,
                                               const CBasisFunction &a_function,
                                               const CBasisFunction &b_function,
                                               const CBasisFunction &c_function,
                                               const CSimdMatrix    &coordinates,
                                               const CSimdMatrix    &c_coordinates,
                                               CSimdMatrix          &buffer,
                                               const double          omega,
                                               const double          threshold) -> void
{
    if (npairs > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_rs_ifl_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 407720, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 3094 * natoms * npairs, 0.0);

        return;
    }

    const auto pi = mathconst::pi_value();

    // NOTE: a row of the values spans every atom pair of every atom on the ket
    // side, so a kernel handed the block of one atom steps by this to reach the
    // next component -- which is what lets it be the kernel a two-center form
    // uses, unchanged.

    const auto nvalues = natoms * npairs;

    simdfunc::compute_pair_exponents(a_function, b_function, coordinates, nmax);

    for (size_t n = 0; n < natoms; n++)
    {
        simdfunc::prepare_buffer(buffer, 407720, 350358, 19401, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 17,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 25, 3, 17,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 44, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 47, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 50, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 53, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 56, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 59, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 62, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 65, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 68, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 71, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 74, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 77, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 80, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 83, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 86, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 89, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 92, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 95, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 98, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 101, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 104, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 107, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 110, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 113, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 116, 0, 3, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 119, 0, 3, 34, 35,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 122, 0, 3, 35, 36,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 125, 0, 3, 36, 37,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 128, 0, 3, 37, 38,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 131, 0, 3, 38, 39,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 134, 0, 3, 39, 40,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 137, 0, 3, 40, 41,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 140, 0, 3, 41, 42,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 143, 0, 3, 42, 43,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 146, 0, 3, 7, 8,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 152, 0, 3, 8, 9,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 158, 0, 3, 9, 10,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 164, 0, 3, 10, 11,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 170, 0, 3, 11, 12,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 176, 0, 3, 12, 13,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 182, 0, 3, 13, 14,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 188, 0, 3, 14, 15,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 194, 0, 3, 15, 16,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 200, 0, 3, 16, 17,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 206, 0, 3, 17, 18,
                                                                       74, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 212, 0, 3, 18, 19,
                                                                       77, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 218, 0, 3, 19, 20,
                                                                       80, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 224, 0, 3, 20, 21,
                                                                       83, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 230, 0, 3, 21, 22,
                                                                       86, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 236, 0, 3, 22, 23,
                                                                       89, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 242, 0, 3, 26, 27,
                                                                       95, 98, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 248, 0, 3, 27, 28,
                                                                       98, 101, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 254, 0, 3, 28, 29,
                                                                       101, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 260, 0, 3, 29, 30,
                                                                       104, 107, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 266, 0, 3, 30, 31,
                                                                       107, 110, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 272, 0, 3, 31, 32,
                                                                       110, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 278, 0, 3, 32, 33,
                                                                       113, 116, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 284, 0, 3, 33, 34,
                                                                       116, 119, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 290, 0, 3, 34, 35,
                                                                       119, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 296, 0, 3, 35, 36,
                                                                       122, 125, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 302, 0, 3, 36, 37,
                                                                       125, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 308, 0, 3, 37, 38,
                                                                       128, 131, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 314, 0, 3, 38, 39,
                                                                       131, 134, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 320, 0, 3, 39, 40,
                                                                       134, 137, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 326, 0, 3, 40, 41,
                                                                       137, 140, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 332, 0, 3, 41, 42,
                                                                       140, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 338, 0, 3, 44, 47,
                                                                       146, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 348, 0, 3, 47, 50,
                                                                       152, 158, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 50, 53,
                                                                       158, 164, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 368, 0, 3, 53, 56,
                                                                       164, 170, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 378, 0, 3, 56, 59,
                                                                       170, 176, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 388, 0, 3, 59, 62,
                                                                       176, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 398, 0, 3, 62, 65,
                                                                       182, 188, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 408, 0, 3, 65, 68,
                                                                       188, 194, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 418, 0, 3, 68, 71,
                                                                       194, 200, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 71, 74,
                                                                       200, 206, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 438, 0, 3, 74, 77,
                                                                       206, 212, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 448, 0, 3, 77, 80,
                                                                       212, 218, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 458, 0, 3, 80, 83,
                                                                       218, 224, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 468, 0, 3, 83, 86,
                                                                       224, 230, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 478, 0, 3, 86, 89,
                                                                       230, 236, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 488, 0, 3, 95, 98,
                                                                       242, 248, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 498, 0, 3, 98,
                                                                       101, 248, 254, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 508, 0, 3, 101,
                                                                       104, 254, 260, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 518, 0, 3, 104,
                                                                       107, 260, 266, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 528, 0, 3, 107,
                                                                       110, 266, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 538, 0, 3, 110,
                                                                       113, 272, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 548, 0, 3, 113,
                                                                       116, 278, 284, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 558, 0, 3, 116,
                                                                       119, 284, 290, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 568, 0, 3, 119,
                                                                       122, 290, 296, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 578, 0, 3, 122,
                                                                       125, 296, 302, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 588, 0, 3, 125,
                                                                       128, 302, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 598, 0, 3, 128,
                                                                       131, 308, 314, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 608, 0, 3, 131,
                                                                       134, 314, 320, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 618, 0, 3, 134,
                                                                       137, 320, 326, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 628, 0, 3, 137,
                                                                       140, 326, 332, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 638, 0, 3, 146,
                                                                       152, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 653, 0, 3, 152,
                                                                       158, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 668, 0, 3, 158,
                                                                       164, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 683, 0, 3, 164,
                                                                       170, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 698, 0, 3, 170,
                                                                       176, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 713, 0, 3, 176,
                                                                       182, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 728, 0, 3, 182,
                                                                       188, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 743, 0, 3, 188,
                                                                       194, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 758, 0, 3, 194,
                                                                       200, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 773, 0, 3, 200,
                                                                       206, 428, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 788, 0, 3, 206,
                                                                       212, 438, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 803, 0, 3, 212,
                                                                       218, 448, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 818, 0, 3, 218,
                                                                       224, 458, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 833, 0, 3, 224,
                                                                       230, 468, 478, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 848, 0, 3, 242,
                                                                       248, 488, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 863, 0, 3, 248,
                                                                       254, 498, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 878, 0, 3, 254,
                                                                       260, 508, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 893, 0, 3, 260,
                                                                       266, 518, 528, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 908, 0, 3, 266,
                                                                       272, 528, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 923, 0, 3, 272,
                                                                       278, 538, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 938, 0, 3, 278,
                                                                       284, 548, 558, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 953, 0, 3, 284,
                                                                       290, 558, 568, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 968, 0, 3, 290,
                                                                       296, 568, 578, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 983, 0, 3, 296,
                                                                       302, 578, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 998, 0, 3, 302,
                                                                       308, 588, 598, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1013, 0, 3, 308,
                                                                       314, 598, 608, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1028, 0, 3, 314,
                                                                       320, 608, 618, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1043, 0, 3, 320,
                                                                       326, 618, 628, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1058, 0, 3, 338,
                                                                       348, 638, 653, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1079, 0, 3, 348,
                                                                       358, 653, 668, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1100, 0, 3, 358,
                                                                       368, 668, 683, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1121, 0, 3, 368,
                                                                       378, 683, 698, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1142, 0, 3, 378,
                                                                       388, 698, 713, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1163, 0, 3, 388,
                                                                       398, 713, 728, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1184, 0, 3, 398,
                                                                       408, 728, 743, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1205, 0, 3, 408,
                                                                       418, 743, 758, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1226, 0, 3, 418,
                                                                       428, 758, 773, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1247, 0, 3, 428,
                                                                       438, 773, 788, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 438,
                                                                       448, 788, 803, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1289, 0, 3, 448,
                                                                       458, 803, 818, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1310, 0, 3, 458,
                                                                       468, 818, 833, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1331, 0, 3, 488,
                                                                       498, 848, 863, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1352, 0, 3, 498,
                                                                       508, 863, 878, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1373, 0, 3, 508,
                                                                       518, 878, 893, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1394, 0, 3, 518,
                                                                       528, 893, 908, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1415, 0, 3, 528,
                                                                       538, 908, 923, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1436, 0, 3, 538,
                                                                       548, 923, 938, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1457, 0, 3, 548,
                                                                       558, 938, 953, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1478, 0, 3, 558,
                                                                       568, 953, 968, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1499, 0, 3, 568,
                                                                       578, 968, 983, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 578,
                                                                       588, 983, 998, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1541, 0, 3, 588,
                                                                       598, 998, 1013, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1562, 0, 3, 598,
                                                                       608, 1013, 1028, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1583, 0, 3, 608,
                                                                       618, 1028, 1043, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1604, 0, 3, 638,
                                                                       653, 1058, 1079, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1632, 0, 3, 653,
                                                                       668, 1079, 1100, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1660, 0, 3, 668,
                                                                       683, 1100, 1121, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 683,
                                                                       698, 1121, 1142, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1716, 0, 3, 698,
                                                                       713, 1142, 1163, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1744, 0, 3, 713,
                                                                       728, 1163, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1772, 0, 3, 728,
                                                                       743, 1184, 1205, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1800, 0, 3, 743,
                                                                       758, 1205, 1226, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1828, 0, 3, 758,
                                                                       773, 1226, 1247, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1856, 0, 3, 773,
                                                                       788, 1247, 1268, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1884, 0, 3, 788,
                                                                       803, 1268, 1289, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1912, 0, 3, 803,
                                                                       818, 1289, 1310, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1940, 0, 3, 848,
                                                                       863, 1331, 1352, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1968, 0, 3, 863,
                                                                       878, 1352, 1373, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1996, 0, 3, 878,
                                                                       893, 1373, 1394, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2024, 0, 3, 893,
                                                                       908, 1394, 1415, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2052, 0, 3, 908,
                                                                       923, 1415, 1436, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2080, 0, 3, 923,
                                                                       938, 1436, 1457, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2108, 0, 3, 938,
                                                                       953, 1457, 1478, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2136, 0, 3, 953,
                                                                       968, 1478, 1499, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2164, 0, 3, 968,
                                                                       983, 1499, 1520, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2192, 0, 3, 983,
                                                                       998, 1520, 1541, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2220, 0, 3, 998,
                                                                       1013, 1541, 1562, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2248, 0, 3, 1013,
                                                                       1028, 1562, 1583, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2276, 0, 3, 1058,
                                                                       1079, 1604, 1632, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2312, 0, 3, 1079,
                                                                       1100, 1632, 1660, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2348, 0, 3, 1100,
                                                                       1121, 1660, 1688, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2384, 0, 3, 1121,
                                                                       1142, 1688, 1716, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2420, 0, 3, 1142,
                                                                       1163, 1716, 1744, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2456, 0, 3, 1163,
                                                                       1184, 1744, 1772, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2492, 0, 3, 1184,
                                                                       1205, 1772, 1800, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2528, 0, 3, 1205,
                                                                       1226, 1800, 1828, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2564, 0, 3, 1226,
                                                                       1247, 1828, 1856, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2600, 0, 3, 1247,
                                                                       1268, 1856, 1884, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2636, 0, 3, 1268,
                                                                       1289, 1884, 1912, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2672, 0, 3, 1331,
                                                                       1352, 1940, 1968, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2708, 0, 3, 1352,
                                                                       1373, 1968, 1996, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2744, 0, 3, 1373,
                                                                       1394, 1996, 2024, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2780, 0, 3, 1394,
                                                                       1415, 2024, 2052, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2816, 0, 3, 1415,
                                                                       1436, 2052, 2080, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2852, 0, 3, 1436,
                                                                       1457, 2080, 2108, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2888, 0, 3, 1457,
                                                                       1478, 2108, 2136, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2924, 0, 3, 1478,
                                                                       1499, 2136, 2164, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2960, 0, 3, 1499,
                                                                       1520, 2164, 2192, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2996, 0, 3, 1520,
                                                                       1541, 2192, 2220, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3032, 0, 3, 1541,
                                                                       1562, 2220, 2248, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3068, 0, 3, 1604,
                                                                       1632, 2276, 2312, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3113, 0, 3, 1632,
                                                                       1660, 2312, 2348, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3158, 0, 3, 1660,
                                                                       1688, 2348, 2384, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3203, 0, 3, 1688,
                                                                       1716, 2384, 2420, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3248, 0, 3, 1716,
                                                                       1744, 2420, 2456, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3293, 0, 3, 1744,
                                                                       1772, 2456, 2492, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3338, 0, 3, 1772,
                                                                       1800, 2492, 2528, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3383, 0, 3, 1800,
                                                                       1828, 2528, 2564, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3428, 0, 3, 1828,
                                                                       1856, 2564, 2600, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3473, 0, 3, 1856,
                                                                       1884, 2600, 2636, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3518, 0, 3, 1940,
                                                                       1968, 2672, 2708, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3563, 0, 3, 1968,
                                                                       1996, 2708, 2744, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3608, 0, 3, 1996,
                                                                       2024, 2744, 2780, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3653, 0, 3, 2024,
                                                                       2052, 2780, 2816, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3698, 0, 3, 2052,
                                                                       2080, 2816, 2852, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3743, 0, 3, 2080,
                                                                       2108, 2852, 2888, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3788, 0, 3, 2108,
                                                                       2136, 2888, 2924, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3833, 0, 3, 2136,
                                                                       2164, 2924, 2960, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3878, 0, 3, 2164,
                                                                       2192, 2960, 2996, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3923, 0, 3, 2192,
                                                                       2220, 2996, 3032, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3968, 0, 3, 2276,
                                                                       2312, 3068, 3113, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4023, 0, 3, 2312,
                                                                       2348, 3113, 3158, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4078, 0, 3, 2348,
                                                                       2384, 3158, 3203, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4133, 0, 3, 2384,
                                                                       2420, 3203, 3248, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4188, 0, 3, 2420,
                                                                       2456, 3248, 3293, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4243, 0, 3, 2456,
                                                                       2492, 3293, 3338, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4298, 0, 3, 2492,
                                                                       2528, 3338, 3383, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4353, 0, 3, 2528,
                                                                       2564, 3383, 3428, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4408, 0, 3, 2564,
                                                                       2600, 3428, 3473, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4463, 0, 3, 2672,
                                                                       2708, 3518, 3563, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4518, 0, 3, 2708,
                                                                       2744, 3563, 3608, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4573, 0, 3, 2744,
                                                                       2780, 3608, 3653, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4628, 0, 3, 2780,
                                                                       2816, 3653, 3698, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4683, 0, 3, 2816,
                                                                       2852, 3698, 3743, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4738, 0, 3, 2852,
                                                                       2888, 3743, 3788, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4793, 0, 3, 2888,
                                                                       2924, 3788, 3833, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4848, 0, 3, 2924,
                                                                       2960, 3833, 3878, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4903, 0, 3, 2960,
                                                                       2996, 3878, 3923, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4958, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4961, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4964, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4967, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4970, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4973, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4976, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4979, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4982, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4985, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4988, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4991, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4994, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4997, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5000, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5003, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5006, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5009, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5012, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5015, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5018, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5021, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5024, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5027, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5030, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5033, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5036, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5039, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5042, 3, 40,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5045, 3, 41,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5048, 3, 42,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5051, 3, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5054, 3, 9, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5063, 3, 10, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5072, 3, 11, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5081, 3, 12, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5090, 3, 13, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5099, 3, 14, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5108, 3, 15, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5117, 3, 16, 71,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5126, 3, 17, 74,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5135, 3, 18, 77,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5144, 3, 19, 80,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5153, 3, 20, 83,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5162, 3, 21, 86,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5171, 3, 22, 89,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5180, 3, 23, 92,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5189, 3, 28, 101,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5198, 3, 29, 104,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5207, 3, 30, 107,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5216, 3, 31, 110,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5225, 3, 32, 113,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5234, 3, 33, 116,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5243, 3, 34, 119,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5252, 3, 35, 122,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5261, 3, 36, 125,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5270, 3, 37, 128,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5279, 3, 38, 131,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5288, 3, 39, 134,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5297, 3, 40, 137,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5306, 3, 41, 140,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5315, 3, 42, 143,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5324, 3, 50, 158,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5342, 3, 53, 164,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5360, 3, 56, 170,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5378, 3, 59, 176,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5396, 3, 62, 182,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5414, 3, 65, 188,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5432, 3, 68, 194,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5450, 3, 71, 200,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5468, 3, 74, 206,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5486, 3, 77, 212,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5504, 3, 80, 218,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5522, 3, 83, 224,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5540, 3, 86, 230,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5558, 3, 89, 236,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5576, 3, 101, 254,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5594, 3, 104, 260,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5612, 3, 107, 266,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5630, 3, 110, 272,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5648, 3, 113, 278,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5666, 3, 116, 284,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5684, 3, 119, 290,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5702, 3, 122, 296,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5720, 3, 125, 302,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5738, 3, 128, 308,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5756, 3, 131, 314,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5774, 3, 134, 320,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5792, 3, 137, 326,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5810, 3, 140, 332,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5828, 3, 158, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5858, 3, 164, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5888, 3, 170, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5918, 3, 176, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5948, 3, 182, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5978, 3, 188, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6008, 3, 194, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6038, 3, 200, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6068, 3, 206, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6098, 3, 212, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6128, 3, 218, 458,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6158, 3, 224, 468,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6188, 3, 230, 478,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6218, 3, 254, 508,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6248, 3, 260, 518,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6278, 3, 266, 528,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6308, 3, 272, 538,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6338, 3, 278, 548,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6368, 3, 284, 558,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6398, 3, 290, 568,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6428, 3, 296, 578,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6458, 3, 302, 588,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6488, 3, 308, 598,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6518, 3, 314, 608,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6548, 3, 320, 618,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6578, 3, 326, 628,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6608, 3, 358, 668,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6653, 3, 368, 683,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6698, 3, 378, 698,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6743, 3, 388, 713,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6788, 3, 398, 728,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6833, 3, 408, 743,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6878, 3, 418, 758,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6923, 3, 428, 773,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6968, 3, 438, 788,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7013, 3, 448, 803,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7058, 3, 458, 818,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7103, 3, 468, 833,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7148, 3, 508, 878,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7193, 3, 518, 893,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7238, 3, 528, 908,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7283, 3, 538, 923,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7328, 3, 548, 938,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7373, 3, 558, 953,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7418, 3, 568, 968,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7463, 3, 578, 983,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7508, 3, 588, 998,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7553, 3, 598,
                                                                       1013, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7598, 3, 608,
                                                                       1028, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7643, 3, 618,
                                                                       1043, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7688, 3, 668,
                                                                       1100, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7751, 3, 683,
                                                                       1121, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7814, 3, 698,
                                                                       1142, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7877, 3, 713,
                                                                       1163, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7940, 3, 728,
                                                                       1184, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8003, 3, 743,
                                                                       1205, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8066, 3, 758,
                                                                       1226, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8129, 3, 773,
                                                                       1247, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8192, 3, 788,
                                                                       1268, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8255, 3, 803,
                                                                       1289, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8318, 3, 818,
                                                                       1310, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8381, 3, 878,
                                                                       1373, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8444, 3, 893,
                                                                       1394, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8507, 3, 908,
                                                                       1415, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8570, 3, 923,
                                                                       1436, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8633, 3, 938,
                                                                       1457, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8696, 3, 953,
                                                                       1478, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8759, 3, 968,
                                                                       1499, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8822, 3, 983,
                                                                       1520, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8885, 3, 998,
                                                                       1541, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8948, 3, 1013,
                                                                       1562, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9011, 3, 1028,
                                                                       1583, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9074, 3, 1100,
                                                                       1660, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9158, 3, 1121,
                                                                       1688, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9242, 3, 1142,
                                                                       1716, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9326, 3, 1163,
                                                                       1744, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9410, 3, 1184,
                                                                       1772, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9494, 3, 1205,
                                                                       1800, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9578, 3, 1226,
                                                                       1828, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9662, 3, 1247,
                                                                       1856, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9746, 3, 1268,
                                                                       1884, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9830, 3, 1289,
                                                                       1912, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9914, 3, 1373,
                                                                       1996, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9998, 3, 1394,
                                                                       2024, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10082, 3, 1415,
                                                                       2052, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10166, 3, 1436,
                                                                       2080, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10250, 3, 1457,
                                                                       2108, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10334, 3, 1478,
                                                                       2136, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10418, 3, 1499,
                                                                       2164, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10502, 3, 1520,
                                                                       2192, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10586, 3, 1541,
                                                                       2220, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10670, 3, 1562,
                                                                       2248, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10754, 3, 1660,
                                                                       2348, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10862, 3, 1688,
                                                                       2384, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10970, 3, 1716,
                                                                       2420, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11078, 3, 1744,
                                                                       2456, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11186, 3, 1772,
                                                                       2492, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11294, 3, 1800,
                                                                       2528, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11402, 3, 1828,
                                                                       2564, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11510, 3, 1856,
                                                                       2600, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11618, 3, 1884,
                                                                       2636, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11726, 3, 1996,
                                                                       2744, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11834, 3, 2024,
                                                                       2780, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11942, 3, 2052,
                                                                       2816, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12050, 3, 2080,
                                                                       2852, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12158, 3, 2108,
                                                                       2888, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12266, 3, 2136,
                                                                       2924, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12374, 3, 2164,
                                                                       2960, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12482, 3, 2192,
                                                                       2996, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12590, 3, 2220,
                                                                       3032, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12698, 3, 2348,
                                                                       3158, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12833, 3, 2384,
                                                                       3203, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12968, 3, 2420,
                                                                       3248, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13103, 3, 2456,
                                                                       3293, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13238, 3, 2492,
                                                                       3338, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13373, 3, 2528,
                                                                       3383, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13508, 3, 2564,
                                                                       3428, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13643, 3, 2600,
                                                                       3473, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13778, 3, 2744,
                                                                       3608, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13913, 3, 2780,
                                                                       3653, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14048, 3, 2816,
                                                                       3698, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14183, 3, 2852,
                                                                       3743, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14318, 3, 2888,
                                                                       3788, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14453, 3, 2924,
                                                                       3833, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14588, 3, 2960,
                                                                       3878, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14723, 3, 2996,
                                                                       3923, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14858, 3, 3158,
                                                                       4078, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15023, 3, 3203,
                                                                       4133, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15188, 3, 3248,
                                                                       4188, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15353, 3, 3293,
                                                                       4243, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15518, 3, 3338,
                                                                       4298, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15683, 3, 3383,
                                                                       4353, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15848, 3, 3428,
                                                                       4408, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16013, 3, 3608,
                                                                       4573, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16178, 3, 3653,
                                                                       4628, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16343, 3, 3698,
                                                                       4683, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16508, 3, 3743,
                                                                       4738, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16673, 3, 3788,
                                                                       4793, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16838, 3, 3833,
                                                                       4848, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17003, 3, 3878,
                                                                       4903, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17168, 3, 7, 8,
                                                                       4958, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17174, 3, 8, 9,
                                                                       4961, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17180, 3, 9, 10,
                                                                       4964, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17186, 3, 10, 11,
                                                                       4967, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17192, 3, 11, 12,
                                                                       4970, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17198, 3, 12, 13,
                                                                       4973, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17204, 3, 13, 14,
                                                                       4976, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17210, 3, 14, 15,
                                                                       4979, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17216, 3, 15, 16,
                                                                       4982, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17222, 3, 16, 17,
                                                                       4985, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17228, 3, 17, 18,
                                                                       4988, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17234, 3, 18, 19,
                                                                       4991, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17240, 3, 19, 20,
                                                                       4994, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17246, 3, 20, 21,
                                                                       4997, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17252, 3, 21, 22,
                                                                       5000, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17258, 3, 22, 23,
                                                                       5003, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17264, 3, 26, 27,
                                                                       5006, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17270, 3, 27, 28,
                                                                       5009, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17276, 3, 28, 29,
                                                                       5012, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17282, 3, 29, 30,
                                                                       5015, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17288, 3, 30, 31,
                                                                       5018, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17294, 3, 31, 32,
                                                                       5021, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17300, 3, 32, 33,
                                                                       5024, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17306, 3, 33, 34,
                                                                       5027, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17312, 3, 34, 35,
                                                                       5030, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17318, 3, 35, 36,
                                                                       5033, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17324, 3, 36, 37,
                                                                       5036, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17330, 3, 37, 38,
                                                                       5039, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17336, 3, 38, 39,
                                                                       5042, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17342, 3, 39, 40,
                                                                       5045, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17348, 3, 40, 41,
                                                                       5048, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17354, 3, 41, 42,
                                                                       5051, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17360, 0, 3,
                                                                       17168, 4958, 17174, 5054,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17378, 0, 3,
                                                                       17174, 4961, 17180, 5063,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17396, 0, 3,
                                                                       17180, 4964, 17186, 5072,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17414, 0, 3,
                                                                       17186, 4967, 17192, 5081,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17432, 0, 3,
                                                                       17192, 4970, 17198, 5090,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17450, 0, 3,
                                                                       17198, 4973, 17204, 5099,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17468, 0, 3,
                                                                       17204, 4976, 17210, 5108,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17486, 0, 3,
                                                                       17210, 4979, 17216, 5117,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17504, 0, 3,
                                                                       17216, 4982, 17222, 5126,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17522, 0, 3,
                                                                       17222, 4985, 17228, 5135,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17540, 0, 3,
                                                                       17228, 4988, 17234, 5144,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17558, 0, 3,
                                                                       17234, 4991, 17240, 5153,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17576, 0, 3,
                                                                       17240, 4994, 17246, 5162,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17594, 0, 3,
                                                                       17246, 4997, 17252, 5171,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17612, 0, 3,
                                                                       17252, 5000, 17258, 5180,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17630, 0, 3,
                                                                       17264, 5006, 17270, 5189,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17648, 0, 3,
                                                                       17270, 5009, 17276, 5198,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17666, 0, 3,
                                                                       17276, 5012, 17282, 5207,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17684, 0, 3,
                                                                       17282, 5015, 17288, 5216,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17702, 0, 3,
                                                                       17288, 5018, 17294, 5225,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17720, 0, 3,
                                                                       17294, 5021, 17300, 5234,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17738, 0, 3,
                                                                       17300, 5024, 17306, 5243,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17756, 0, 3,
                                                                       17306, 5027, 17312, 5252,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17774, 0, 3,
                                                                       17312, 5030, 17318, 5261,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17792, 0, 3,
                                                                       17318, 5033, 17324, 5270,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17810, 0, 3,
                                                                       17324, 5036, 17330, 5279,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17828, 0, 3,
                                                                       17330, 5039, 17336, 5288,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17846, 0, 3,
                                                                       17336, 5042, 17342, 5297,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17864, 0, 3,
                                                                       17342, 5045, 17348, 5306,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17882, 0, 3,
                                                                       17348, 5048, 17354, 5315,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17900, 0, 3,
                                                                       17360, 5054, 17378, 146,
                                                                       152, 5324, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17936, 0, 3,
                                                                       17378, 5063, 17396, 152,
                                                                       158, 5342, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17972, 0, 3,
                                                                       17396, 5072, 17414, 158,
                                                                       164, 5360, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18008, 0, 3,
                                                                       17414, 5081, 17432, 164,
                                                                       170, 5378, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18044, 0, 3,
                                                                       17432, 5090, 17450, 170,
                                                                       176, 5396, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18080, 0, 3,
                                                                       17450, 5099, 17468, 176,
                                                                       182, 5414, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18116, 0, 3,
                                                                       17468, 5108, 17486, 182,
                                                                       188, 5432, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18152, 0, 3,
                                                                       17486, 5117, 17504, 188,
                                                                       194, 5450, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18188, 0, 3,
                                                                       17504, 5126, 17522, 194,
                                                                       200, 5468, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18224, 0, 3,
                                                                       17522, 5135, 17540, 200,
                                                                       206, 5486, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18260, 0, 3,
                                                                       17540, 5144, 17558, 206,
                                                                       212, 5504, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18296, 0, 3,
                                                                       17558, 5153, 17576, 212,
                                                                       218, 5522, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18332, 0, 3,
                                                                       17576, 5162, 17594, 218,
                                                                       224, 5540, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18368, 0, 3,
                                                                       17594, 5171, 17612, 224,
                                                                       230, 5558, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18404, 0, 3,
                                                                       17630, 5189, 17648, 242,
                                                                       248, 5576, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18440, 0, 3,
                                                                       17648, 5198, 17666, 248,
                                                                       254, 5594, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18476, 0, 3,
                                                                       17666, 5207, 17684, 254,
                                                                       260, 5612, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18512, 0, 3,
                                                                       17684, 5216, 17702, 260,
                                                                       266, 5630, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18548, 0, 3,
                                                                       17702, 5225, 17720, 266,
                                                                       272, 5648, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18584, 0, 3,
                                                                       17720, 5234, 17738, 272,
                                                                       278, 5666, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18620, 0, 3,
                                                                       17738, 5243, 17756, 278,
                                                                       284, 5684, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18656, 0, 3,
                                                                       17756, 5252, 17774, 284,
                                                                       290, 5702, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18692, 0, 3,
                                                                       17774, 5261, 17792, 290,
                                                                       296, 5720, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18728, 0, 3,
                                                                       17792, 5270, 17810, 296,
                                                                       302, 5738, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18764, 0, 3,
                                                                       17810, 5279, 17828, 302,
                                                                       308, 5756, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18800, 0, 3,
                                                                       17828, 5288, 17846, 308,
                                                                       314, 5774, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18836, 0, 3,
                                                                       17846, 5297, 17864, 314,
                                                                       320, 5792, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18872, 0, 3,
                                                                       17864, 5306, 17882, 320,
                                                                       326, 5810, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18908, 0, 3,
                                                                       17900, 5324, 17936, 338,
                                                                       348, 5828, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18968, 0, 3,
                                                                       17936, 5342, 17972, 348,
                                                                       358, 5858, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19028, 0, 3,
                                                                       17972, 5360, 18008, 358,
                                                                       368, 5888, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19088, 0, 3,
                                                                       18008, 5378, 18044, 368,
                                                                       378, 5918, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19148, 0, 3,
                                                                       18044, 5396, 18080, 378,
                                                                       388, 5948, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19208, 0, 3,
                                                                       18080, 5414, 18116, 388,
                                                                       398, 5978, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19268, 0, 3,
                                                                       18116, 5432, 18152, 398,
                                                                       408, 6008, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19328, 0, 3,
                                                                       18152, 5450, 18188, 408,
                                                                       418, 6038, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19388, 0, 3,
                                                                       18188, 5468, 18224, 418,
                                                                       428, 6068, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19448, 0, 3,
                                                                       18224, 5486, 18260, 428,
                                                                       438, 6098, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19508, 0, 3,
                                                                       18260, 5504, 18296, 438,
                                                                       448, 6128, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19568, 0, 3,
                                                                       18296, 5522, 18332, 448,
                                                                       458, 6158, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19628, 0, 3,
                                                                       18332, 5540, 18368, 458,
                                                                       468, 6188, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19688, 0, 3,
                                                                       18404, 5576, 18440, 488,
                                                                       498, 6218, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19748, 0, 3,
                                                                       18440, 5594, 18476, 498,
                                                                       508, 6248, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19808, 0, 3,
                                                                       18476, 5612, 18512, 508,
                                                                       518, 6278, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19868, 0, 3,
                                                                       18512, 5630, 18548, 518,
                                                                       528, 6308, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19928, 0, 3,
                                                                       18548, 5648, 18584, 528,
                                                                       538, 6338, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19988, 0, 3,
                                                                       18584, 5666, 18620, 538,
                                                                       548, 6368, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20048, 0, 3,
                                                                       18620, 5684, 18656, 548,
                                                                       558, 6398, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20108, 0, 3,
                                                                       18656, 5702, 18692, 558,
                                                                       568, 6428, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20168, 0, 3,
                                                                       18692, 5720, 18728, 568,
                                                                       578, 6458, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20228, 0, 3,
                                                                       18728, 5738, 18764, 578,
                                                                       588, 6488, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20288, 0, 3,
                                                                       18764, 5756, 18800, 588,
                                                                       598, 6518, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20348, 0, 3,
                                                                       18800, 5774, 18836, 598,
                                                                       608, 6548, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20408, 0, 3,
                                                                       18836, 5792, 18872, 608,
                                                                       618, 6578, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20468, 0, 3,
                                                                       18908, 5828, 18968, 638,
                                                                       653, 6608, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20558, 0, 3,
                                                                       18968, 5858, 19028, 653,
                                                                       668, 6653, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20648, 0, 3,
                                                                       19028, 5888, 19088, 668,
                                                                       683, 6698, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20738, 0, 3,
                                                                       19088, 5918, 19148, 683,
                                                                       698, 6743, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20828, 0, 3,
                                                                       19148, 5948, 19208, 698,
                                                                       713, 6788, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20918, 0, 3,
                                                                       19208, 5978, 19268, 713,
                                                                       728, 6833, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21008, 0, 3,
                                                                       19268, 6008, 19328, 728,
                                                                       743, 6878, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21098, 0, 3,
                                                                       19328, 6038, 19388, 743,
                                                                       758, 6923, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21188, 0, 3,
                                                                       19388, 6068, 19448, 758,
                                                                       773, 6968, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21278, 0, 3,
                                                                       19448, 6098, 19508, 773,
                                                                       788, 7013, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21368, 0, 3,
                                                                       19508, 6128, 19568, 788,
                                                                       803, 7058, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21458, 0, 3,
                                                                       19568, 6158, 19628, 803,
                                                                       818, 7103, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21548, 0, 3,
                                                                       19688, 6218, 19748, 848,
                                                                       863, 7148, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21638, 0, 3,
                                                                       19748, 6248, 19808, 863,
                                                                       878, 7193, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21728, 0, 3,
                                                                       19808, 6278, 19868, 878,
                                                                       893, 7238, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21818, 0, 3,
                                                                       19868, 6308, 19928, 893,
                                                                       908, 7283, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21908, 0, 3,
                                                                       19928, 6338, 19988, 908,
                                                                       923, 7328, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21998, 0, 3,
                                                                       19988, 6368, 20048, 923,
                                                                       938, 7373, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22088, 0, 3,
                                                                       20048, 6398, 20108, 938,
                                                                       953, 7418, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22178, 0, 3,
                                                                       20108, 6428, 20168, 953,
                                                                       968, 7463, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22268, 0, 3,
                                                                       20168, 6458, 20228, 968,
                                                                       983, 7508, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22358, 0, 3,
                                                                       20228, 6488, 20288, 983,
                                                                       998, 7553, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22448, 0, 3,
                                                                       20288, 6518, 20348, 998,
                                                                       1013, 7598, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22538, 0, 3,
                                                                       20348, 6548, 20408, 1013,
                                                                       1028, 7643, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22628, 0, 3,
                                                                       20468, 6608, 20558, 1058,
                                                                       1079, 7688, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22754, 0, 3,
                                                                       20558, 6653, 20648, 1079,
                                                                       1100, 7751, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22880, 0, 3,
                                                                       20648, 6698, 20738, 1100,
                                                                       1121, 7814, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23006, 0, 3,
                                                                       20738, 6743, 20828, 1121,
                                                                       1142, 7877, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23132, 0, 3,
                                                                       20828, 6788, 20918, 1142,
                                                                       1163, 7940, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23258, 0, 3,
                                                                       20918, 6833, 21008, 1163,
                                                                       1184, 8003, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23384, 0, 3,
                                                                       21008, 6878, 21098, 1184,
                                                                       1205, 8066, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23510, 0, 3,
                                                                       21098, 6923, 21188, 1205,
                                                                       1226, 8129, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23636, 0, 3,
                                                                       21188, 6968, 21278, 1226,
                                                                       1247, 8192, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23762, 0, 3,
                                                                       21278, 7013, 21368, 1247,
                                                                       1268, 8255, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23888, 0, 3,
                                                                       21368, 7058, 21458, 1268,
                                                                       1289, 8318, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24014, 0, 3,
                                                                       21548, 7148, 21638, 1331,
                                                                       1352, 8381, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24140, 0, 3,
                                                                       21638, 7193, 21728, 1352,
                                                                       1373, 8444, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24266, 0, 3,
                                                                       21728, 7238, 21818, 1373,
                                                                       1394, 8507, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24392, 0, 3,
                                                                       21818, 7283, 21908, 1394,
                                                                       1415, 8570, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24518, 0, 3,
                                                                       21908, 7328, 21998, 1415,
                                                                       1436, 8633, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24644, 0, 3,
                                                                       21998, 7373, 22088, 1436,
                                                                       1457, 8696, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24770, 0, 3,
                                                                       22088, 7418, 22178, 1457,
                                                                       1478, 8759, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24896, 0, 3,
                                                                       22178, 7463, 22268, 1478,
                                                                       1499, 8822, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 25022, 0, 3,
                                                                       22268, 7508, 22358, 1499,
                                                                       1520, 8885, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 25148, 0, 3,
                                                                       22358, 7553, 22448, 1520,
                                                                       1541, 8948, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 25274, 0, 3,
                                                                       22448, 7598, 22538, 1541,
                                                                       1562, 9011, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25400, 0, 3,
                                                                       22628, 7688, 22754, 1604,
                                                                       1632, 9074, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25568, 0, 3,
                                                                       22754, 7751, 22880, 1632,
                                                                       1660, 9158, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25736, 0, 3,
                                                                       22880, 7814, 23006, 1660,
                                                                       1688, 9242, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25904, 0, 3,
                                                                       23006, 7877, 23132, 1688,
                                                                       1716, 9326, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26072, 0, 3,
                                                                       23132, 7940, 23258, 1716,
                                                                       1744, 9410, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26240, 0, 3,
                                                                       23258, 8003, 23384, 1744,
                                                                       1772, 9494, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26408, 0, 3,
                                                                       23384, 8066, 23510, 1772,
                                                                       1800, 9578, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26576, 0, 3,
                                                                       23510, 8129, 23636, 1800,
                                                                       1828, 9662, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26744, 0, 3,
                                                                       23636, 8192, 23762, 1828,
                                                                       1856, 9746, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26912, 0, 3,
                                                                       23762, 8255, 23888, 1856,
                                                                       1884, 9830, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 27080, 0, 3,
                                                                       24014, 8381, 24140, 1940,
                                                                       1968, 9914, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 27248, 0, 3,
                                                                       24140, 8444, 24266, 1968,
                                                                       1996, 9998, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 27416, 0, 3,
                                                                       24266, 8507, 24392, 1996,
                                                                       2024, 10082, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 27584, 0, 3,
                                                                       24392, 8570, 24518, 2024,
                                                                       2052, 10166, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 27752, 0, 3,
                                                                       24518, 8633, 24644, 2052,
                                                                       2080, 10250, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 27920, 0, 3,
                                                                       24644, 8696, 24770, 2080,
                                                                       2108, 10334, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 28088, 0, 3,
                                                                       24770, 8759, 24896, 2108,
                                                                       2136, 10418, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 28256, 0, 3,
                                                                       24896, 8822, 25022, 2136,
                                                                       2164, 10502, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 28424, 0, 3,
                                                                       25022, 8885, 25148, 2164,
                                                                       2192, 10586, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 28592, 0, 3,
                                                                       25148, 8948, 25274, 2192,
                                                                       2220, 10670, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28760, 0, 3,
                                                                       25400, 9074, 25568, 2276,
                                                                       2312, 10754, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28976, 0, 3,
                                                                       25568, 9158, 25736, 2312,
                                                                       2348, 10862, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 29192, 0, 3,
                                                                       25736, 9242, 25904, 2348,
                                                                       2384, 10970, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 29408, 0, 3,
                                                                       25904, 9326, 26072, 2384,
                                                                       2420, 11078, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 29624, 0, 3,
                                                                       26072, 9410, 26240, 2420,
                                                                       2456, 11186, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 29840, 0, 3,
                                                                       26240, 9494, 26408, 2456,
                                                                       2492, 11294, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 30056, 0, 3,
                                                                       26408, 9578, 26576, 2492,
                                                                       2528, 11402, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 30272, 0, 3,
                                                                       26576, 9662, 26744, 2528,
                                                                       2564, 11510, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 30488, 0, 3,
                                                                       26744, 9746, 26912, 2564,
                                                                       2600, 11618, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 30704, 0, 3,
                                                                       27080, 9914, 27248, 2672,
                                                                       2708, 11726, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 30920, 0, 3,
                                                                       27248, 9998, 27416, 2708,
                                                                       2744, 11834, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 31136, 0, 3,
                                                                       27416, 10082, 27584, 2744,
                                                                       2780, 11942, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 31352, 0, 3,
                                                                       27584, 10166, 27752, 2780,
                                                                       2816, 12050, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 31568, 0, 3,
                                                                       27752, 10250, 27920, 2816,
                                                                       2852, 12158, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 31784, 0, 3,
                                                                       27920, 10334, 28088, 2852,
                                                                       2888, 12266, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 32000, 0, 3,
                                                                       28088, 10418, 28256, 2888,
                                                                       2924, 12374, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 32216, 0, 3,
                                                                       28256, 10502, 28424, 2924,
                                                                       2960, 12482, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 32432, 0, 3,
                                                                       28424, 10586, 28592, 2960,
                                                                       2996, 12590, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 32648, 0, 3,
                                                                       28760, 10754, 28976, 3068,
                                                                       3113, 12698, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 32918, 0, 3,
                                                                       28976, 10862, 29192, 3113,
                                                                       3158, 12833, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 33188, 0, 3,
                                                                       29192, 10970, 29408, 3158,
                                                                       3203, 12968, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 33458, 0, 3,
                                                                       29408, 11078, 29624, 3203,
                                                                       3248, 13103, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 33728, 0, 3,
                                                                       29624, 11186, 29840, 3248,
                                                                       3293, 13238, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 33998, 0, 3,
                                                                       29840, 11294, 30056, 3293,
                                                                       3338, 13373, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 34268, 0, 3,
                                                                       30056, 11402, 30272, 3338,
                                                                       3383, 13508, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 34538, 0, 3,
                                                                       30272, 11510, 30488, 3383,
                                                                       3428, 13643, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 34808, 0, 3,
                                                                       30704, 11726, 30920, 3518,
                                                                       3563, 13778, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 35078, 0, 3,
                                                                       30920, 11834, 31136, 3563,
                                                                       3608, 13913, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 35348, 0, 3,
                                                                       31136, 11942, 31352, 3608,
                                                                       3653, 14048, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 35618, 0, 3,
                                                                       31352, 12050, 31568, 3653,
                                                                       3698, 14183, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 35888, 0, 3,
                                                                       31568, 12158, 31784, 3698,
                                                                       3743, 14318, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 36158, 0, 3,
                                                                       31784, 12266, 32000, 3743,
                                                                       3788, 14453, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 36428, 0, 3,
                                                                       32000, 12374, 32216, 3788,
                                                                       3833, 14588, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 36698, 0, 3,
                                                                       32216, 12482, 32432, 3833,
                                                                       3878, 14723, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 36968, 0, 3,
                                                                       32648, 12698, 32918, 3968,
                                                                       4023, 14858, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 37298, 0, 3,
                                                                       32918, 12833, 33188, 4023,
                                                                       4078, 15023, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 37628, 0, 3,
                                                                       33188, 12968, 33458, 4078,
                                                                       4133, 15188, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 37958, 0, 3,
                                                                       33458, 13103, 33728, 4133,
                                                                       4188, 15353, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 38288, 0, 3,
                                                                       33728, 13238, 33998, 4188,
                                                                       4243, 15518, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 38618, 0, 3,
                                                                       33998, 13373, 34268, 4243,
                                                                       4298, 15683, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 38948, 0, 3,
                                                                       34268, 13508, 34538, 4298,
                                                                       4353, 15848, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 39278, 0, 3,
                                                                       34808, 13778, 35078, 4463,
                                                                       4518, 16013, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 39608, 0, 3,
                                                                       35078, 13913, 35348, 4518,
                                                                       4573, 16178, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 39938, 0, 3,
                                                                       35348, 14048, 35618, 4573,
                                                                       4628, 16343, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 40268, 0, 3,
                                                                       35618, 14183, 35888, 4628,
                                                                       4683, 16508, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 40598, 0, 3,
                                                                       35888, 14318, 36158, 4683,
                                                                       4738, 16673, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 40928, 0, 3,
                                                                       36158, 14453, 36428, 4738,
                                                                       4793, 16838, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 41258, 0, 3,
                                                                       36428, 14588, 36698, 4793,
                                                                       4848, 17003, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41588, 3, 4958,
                                                                       4961, 17180, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41598, 3, 4961,
                                                                       4964, 17186, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41608, 3, 4964,
                                                                       4967, 17192, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41618, 3, 4967,
                                                                       4970, 17198, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41628, 3, 4970,
                                                                       4973, 17204, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41638, 3, 4973,
                                                                       4976, 17210, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41648, 3, 4976,
                                                                       4979, 17216, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41658, 3, 4979,
                                                                       4982, 17222, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41668, 3, 4982,
                                                                       4985, 17228, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41678, 3, 4985,
                                                                       4988, 17234, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41688, 3, 4988,
                                                                       4991, 17240, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41698, 3, 4991,
                                                                       4994, 17246, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41708, 3, 4994,
                                                                       4997, 17252, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41718, 3, 4997,
                                                                       5000, 17258, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41728, 3, 5006,
                                                                       5009, 17276, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41738, 3, 5009,
                                                                       5012, 17282, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41748, 3, 5012,
                                                                       5015, 17288, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41758, 3, 5015,
                                                                       5018, 17294, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41768, 3, 5018,
                                                                       5021, 17300, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41778, 3, 5021,
                                                                       5024, 17306, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41788, 3, 5024,
                                                                       5027, 17312, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41798, 3, 5027,
                                                                       5030, 17318, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41808, 3, 5030,
                                                                       5033, 17324, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41818, 3, 5033,
                                                                       5036, 17330, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41828, 3, 5036,
                                                                       5039, 17336, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41838, 3, 5039,
                                                                       5042, 17342, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41848, 3, 5042,
                                                                       5045, 17348, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 41858, 3, 5045,
                                                                       5048, 17354, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 41868, 0, 3,
                                                                       41588, 17180, 41598, 5054,
                                                                       5063, 17396, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 41898, 0, 3,
                                                                       41598, 17186, 41608, 5063,
                                                                       5072, 17414, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 41928, 0, 3,
                                                                       41608, 17192, 41618, 5072,
                                                                       5081, 17432, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 41958, 0, 3,
                                                                       41618, 17198, 41628, 5081,
                                                                       5090, 17450, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 41988, 0, 3,
                                                                       41628, 17204, 41638, 5090,
                                                                       5099, 17468, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42018, 0, 3,
                                                                       41638, 17210, 41648, 5099,
                                                                       5108, 17486, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42048, 0, 3,
                                                                       41648, 17216, 41658, 5108,
                                                                       5117, 17504, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42078, 0, 3,
                                                                       41658, 17222, 41668, 5117,
                                                                       5126, 17522, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42108, 0, 3,
                                                                       41668, 17228, 41678, 5126,
                                                                       5135, 17540, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42138, 0, 3,
                                                                       41678, 17234, 41688, 5135,
                                                                       5144, 17558, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42168, 0, 3,
                                                                       41688, 17240, 41698, 5144,
                                                                       5153, 17576, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42198, 0, 3,
                                                                       41698, 17246, 41708, 5153,
                                                                       5162, 17594, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42228, 0, 3,
                                                                       41708, 17252, 41718, 5162,
                                                                       5171, 17612, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42258, 0, 3,
                                                                       41728, 17276, 41738, 5189,
                                                                       5198, 17666, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42288, 0, 3,
                                                                       41738, 17282, 41748, 5198,
                                                                       5207, 17684, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42318, 0, 3,
                                                                       41748, 17288, 41758, 5207,
                                                                       5216, 17702, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42348, 0, 3,
                                                                       41758, 17294, 41768, 5216,
                                                                       5225, 17720, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42378, 0, 3,
                                                                       41768, 17300, 41778, 5225,
                                                                       5234, 17738, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42408, 0, 3,
                                                                       41778, 17306, 41788, 5234,
                                                                       5243, 17756, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42438, 0, 3,
                                                                       41788, 17312, 41798, 5243,
                                                                       5252, 17774, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42468, 0, 3,
                                                                       41798, 17318, 41808, 5252,
                                                                       5261, 17792, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42498, 0, 3,
                                                                       41808, 17324, 41818, 5261,
                                                                       5270, 17810, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42528, 0, 3,
                                                                       41818, 17330, 41828, 5270,
                                                                       5279, 17828, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42558, 0, 3,
                                                                       41828, 17336, 41838, 5279,
                                                                       5288, 17846, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42588, 0, 3,
                                                                       41838, 17342, 41848, 5288,
                                                                       5297, 17864, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 42618, 0, 3,
                                                                       41848, 17348, 41858, 5297,
                                                                       5306, 17882, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 42648, 0, 3,
                                                                       41868, 17396, 41898, 5324,
                                                                       5342, 17972, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 42708, 0, 3,
                                                                       41898, 17414, 41928, 5342,
                                                                       5360, 18008, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 42768, 0, 3,
                                                                       41928, 17432, 41958, 5360,
                                                                       5378, 18044, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 42828, 0, 3,
                                                                       41958, 17450, 41988, 5378,
                                                                       5396, 18080, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 42888, 0, 3,
                                                                       41988, 17468, 42018, 5396,
                                                                       5414, 18116, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 42948, 0, 3,
                                                                       42018, 17486, 42048, 5414,
                                                                       5432, 18152, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 43008, 0, 3,
                                                                       42048, 17504, 42078, 5432,
                                                                       5450, 18188, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 43068, 0, 3,
                                                                       42078, 17522, 42108, 5450,
                                                                       5468, 18224, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 43128, 0, 3,
                                                                       42108, 17540, 42138, 5468,
                                                                       5486, 18260, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 43188, 0, 3,
                                                                       42138, 17558, 42168, 5486,
                                                                       5504, 18296, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 43248, 0, 3,
                                                                       42168, 17576, 42198, 5504,
                                                                       5522, 18332, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 43308, 0, 3,
                                                                       42198, 17594, 42228, 5522,
                                                                       5540, 18368, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 43368, 0, 3,
                                                                       42258, 17666, 42288, 5576,
                                                                       5594, 18476, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 43428, 0, 3,
                                                                       42288, 17684, 42318, 5594,
                                                                       5612, 18512, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 43488, 0, 3,
                                                                       42318, 17702, 42348, 5612,
                                                                       5630, 18548, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 43548, 0, 3,
                                                                       42348, 17720, 42378, 5630,
                                                                       5648, 18584, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 43608, 0, 3,
                                                                       42378, 17738, 42408, 5648,
                                                                       5666, 18620, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 43668, 0, 3,
                                                                       42408, 17756, 42438, 5666,
                                                                       5684, 18656, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 43728, 0, 3,
                                                                       42438, 17774, 42468, 5684,
                                                                       5702, 18692, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 43788, 0, 3,
                                                                       42468, 17792, 42498, 5702,
                                                                       5720, 18728, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 43848, 0, 3,
                                                                       42498, 17810, 42528, 5720,
                                                                       5738, 18764, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 43908, 0, 3,
                                                                       42528, 17828, 42558, 5738,
                                                                       5756, 18800, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 43968, 0, 3,
                                                                       42558, 17846, 42588, 5756,
                                                                       5774, 18836, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 44028, 0, 3,
                                                                       42588, 17864, 42618, 5774,
                                                                       5792, 18872, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 44088, 0, 3,
                                                                       42648, 17972, 42708, 5828,
                                                                       5858, 19028, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 44188, 0, 3,
                                                                       42708, 18008, 42768, 5858,
                                                                       5888, 19088, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 44288, 0, 3,
                                                                       42768, 18044, 42828, 5888,
                                                                       5918, 19148, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 44388, 0, 3,
                                                                       42828, 18080, 42888, 5918,
                                                                       5948, 19208, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 44488, 0, 3,
                                                                       42888, 18116, 42948, 5948,
                                                                       5978, 19268, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 44588, 0, 3,
                                                                       42948, 18152, 43008, 5978,
                                                                       6008, 19328, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 44688, 0, 3,
                                                                       43008, 18188, 43068, 6008,
                                                                       6038, 19388, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 44788, 0, 3,
                                                                       43068, 18224, 43128, 6038,
                                                                       6068, 19448, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 44888, 0, 3,
                                                                       43128, 18260, 43188, 6068,
                                                                       6098, 19508, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 44988, 0, 3,
                                                                       43188, 18296, 43248, 6098,
                                                                       6128, 19568, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 45088, 0, 3,
                                                                       43248, 18332, 43308, 6128,
                                                                       6158, 19628, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 45188, 0, 3,
                                                                       43368, 18476, 43428, 6218,
                                                                       6248, 19808, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 45288, 0, 3,
                                                                       43428, 18512, 43488, 6248,
                                                                       6278, 19868, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 45388, 0, 3,
                                                                       43488, 18548, 43548, 6278,
                                                                       6308, 19928, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 45488, 0, 3,
                                                                       43548, 18584, 43608, 6308,
                                                                       6338, 19988, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 45588, 0, 3,
                                                                       43608, 18620, 43668, 6338,
                                                                       6368, 20048, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 45688, 0, 3,
                                                                       43668, 18656, 43728, 6368,
                                                                       6398, 20108, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 45788, 0, 3,
                                                                       43728, 18692, 43788, 6398,
                                                                       6428, 20168, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 45888, 0, 3,
                                                                       43788, 18728, 43848, 6428,
                                                                       6458, 20228, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 45988, 0, 3,
                                                                       43848, 18764, 43908, 6458,
                                                                       6488, 20288, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46088, 0, 3,
                                                                       43908, 18800, 43968, 6488,
                                                                       6518, 20348, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46188, 0, 3,
                                                                       43968, 18836, 44028, 6518,
                                                                       6548, 20408, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 46288, 0, 3,
                                                                       44088, 19028, 44188, 6608,
                                                                       6653, 20648, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 46438, 0, 3,
                                                                       44188, 19088, 44288, 6653,
                                                                       6698, 20738, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 46588, 0, 3,
                                                                       44288, 19148, 44388, 6698,
                                                                       6743, 20828, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 46738, 0, 3,
                                                                       44388, 19208, 44488, 6743,
                                                                       6788, 20918, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 46888, 0, 3,
                                                                       44488, 19268, 44588, 6788,
                                                                       6833, 21008, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 47038, 0, 3,
                                                                       44588, 19328, 44688, 6833,
                                                                       6878, 21098, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 47188, 0, 3,
                                                                       44688, 19388, 44788, 6878,
                                                                       6923, 21188, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 47338, 0, 3,
                                                                       44788, 19448, 44888, 6923,
                                                                       6968, 21278, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 47488, 0, 3,
                                                                       44888, 19508, 44988, 6968,
                                                                       7013, 21368, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 47638, 0, 3,
                                                                       44988, 19568, 45088, 7013,
                                                                       7058, 21458, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 47788, 0, 3,
                                                                       45188, 19808, 45288, 7148,
                                                                       7193, 21728, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 47938, 0, 3,
                                                                       45288, 19868, 45388, 7193,
                                                                       7238, 21818, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48088, 0, 3,
                                                                       45388, 19928, 45488, 7238,
                                                                       7283, 21908, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48238, 0, 3,
                                                                       45488, 19988, 45588, 7283,
                                                                       7328, 21998, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48388, 0, 3,
                                                                       45588, 20048, 45688, 7328,
                                                                       7373, 22088, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48538, 0, 3,
                                                                       45688, 20108, 45788, 7373,
                                                                       7418, 22178, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48688, 0, 3,
                                                                       45788, 20168, 45888, 7418,
                                                                       7463, 22268, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48838, 0, 3,
                                                                       45888, 20228, 45988, 7463,
                                                                       7508, 22358, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48988, 0, 3,
                                                                       45988, 20288, 46088, 7508,
                                                                       7553, 22448, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49138, 0, 3,
                                                                       46088, 20348, 46188, 7553,
                                                                       7598, 22538, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 49288, 0, 3,
                                                                       46288, 20648, 46438, 7688,
                                                                       7751, 22880, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 49498, 0, 3,
                                                                       46438, 20738, 46588, 7751,
                                                                       7814, 23006, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 49708, 0, 3,
                                                                       46588, 20828, 46738, 7814,
                                                                       7877, 23132, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 49918, 0, 3,
                                                                       46738, 20918, 46888, 7877,
                                                                       7940, 23258, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 50128, 0, 3,
                                                                       46888, 21008, 47038, 7940,
                                                                       8003, 23384, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 50338, 0, 3,
                                                                       47038, 21098, 47188, 8003,
                                                                       8066, 23510, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 50548, 0, 3,
                                                                       47188, 21188, 47338, 8066,
                                                                       8129, 23636, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 50758, 0, 3,
                                                                       47338, 21278, 47488, 8129,
                                                                       8192, 23762, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 50968, 0, 3,
                                                                       47488, 21368, 47638, 8192,
                                                                       8255, 23888, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51178, 0, 3,
                                                                       47788, 21728, 47938, 8381,
                                                                       8444, 24266, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51388, 0, 3,
                                                                       47938, 21818, 48088, 8444,
                                                                       8507, 24392, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51598, 0, 3,
                                                                       48088, 21908, 48238, 8507,
                                                                       8570, 24518, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51808, 0, 3,
                                                                       48238, 21998, 48388, 8570,
                                                                       8633, 24644, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52018, 0, 3,
                                                                       48388, 22088, 48538, 8633,
                                                                       8696, 24770, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52228, 0, 3,
                                                                       48538, 22178, 48688, 8696,
                                                                       8759, 24896, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52438, 0, 3,
                                                                       48688, 22268, 48838, 8759,
                                                                       8822, 25022, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52648, 0, 3,
                                                                       48838, 22358, 48988, 8822,
                                                                       8885, 25148, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52858, 0, 3,
                                                                       48988, 22448, 49138, 8885,
                                                                       8948, 25274, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 53068, 0, 3,
                                                                       49288, 22880, 49498, 9074,
                                                                       9158, 25736, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 53348, 0, 3,
                                                                       49498, 23006, 49708, 9158,
                                                                       9242, 25904, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 53628, 0, 3,
                                                                       49708, 23132, 49918, 9242,
                                                                       9326, 26072, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 53908, 0, 3,
                                                                       49918, 23258, 50128, 9326,
                                                                       9410, 26240, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54188, 0, 3,
                                                                       50128, 23384, 50338, 9410,
                                                                       9494, 26408, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54468, 0, 3,
                                                                       50338, 23510, 50548, 9494,
                                                                       9578, 26576, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54748, 0, 3,
                                                                       50548, 23636, 50758, 9578,
                                                                       9662, 26744, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 55028, 0, 3,
                                                                       50758, 23762, 50968, 9662,
                                                                       9746, 26912, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 55308, 0, 3,
                                                                       51178, 24266, 51388, 9914,
                                                                       9998, 27416, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 55588, 0, 3,
                                                                       51388, 24392, 51598, 9998,
                                                                       10082, 27584, ncols,
                                                                       gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 55868, 0, 3,
                                                                       51598, 24518, 51808,
                                                                       10082, 10166, 27752,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 56148, 0, 3,
                                                                       51808, 24644, 52018,
                                                                       10166, 10250, 27920,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 56428, 0, 3,
                                                                       52018, 24770, 52228,
                                                                       10250, 10334, 28088,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 56708, 0, 3,
                                                                       52228, 24896, 52438,
                                                                       10334, 10418, 28256,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 56988, 0, 3,
                                                                       52438, 25022, 52648,
                                                                       10418, 10502, 28424,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 57268, 0, 3,
                                                                       52648, 25148, 52858,
                                                                       10502, 10586, 28592,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 57548, 0, 3,
                                                                       53068, 25736, 53348,
                                                                       10754, 10862, 29192,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 57908, 0, 3,
                                                                       53348, 25904, 53628,
                                                                       10862, 10970, 29408,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 58268, 0, 3,
                                                                       53628, 26072, 53908,
                                                                       10970, 11078, 29624,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 58628, 0, 3,
                                                                       53908, 26240, 54188,
                                                                       11078, 11186, 29840,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 58988, 0, 3,
                                                                       54188, 26408, 54468,
                                                                       11186, 11294, 30056,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 59348, 0, 3,
                                                                       54468, 26576, 54748,
                                                                       11294, 11402, 30272,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 59708, 0, 3,
                                                                       54748, 26744, 55028,
                                                                       11402, 11510, 30488,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 60068, 0, 3,
                                                                       55308, 27416, 55588,
                                                                       11726, 11834, 31136,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 60428, 0, 3,
                                                                       55588, 27584, 55868,
                                                                       11834, 11942, 31352,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 60788, 0, 3,
                                                                       55868, 27752, 56148,
                                                                       11942, 12050, 31568,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 61148, 0, 3,
                                                                       56148, 27920, 56428,
                                                                       12050, 12158, 31784,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 61508, 0, 3,
                                                                       56428, 28088, 56708,
                                                                       12158, 12266, 32000,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 61868, 0, 3,
                                                                       56708, 28256, 56988,
                                                                       12266, 12374, 32216,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 62228, 0, 3,
                                                                       56988, 28424, 57268,
                                                                       12374, 12482, 32432,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 62588, 0, 3,
                                                                       57548, 29192, 57908,
                                                                       12698, 12833, 33188,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 63038, 0, 3,
                                                                       57908, 29408, 58268,
                                                                       12833, 12968, 33458,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 63488, 0, 3,
                                                                       58268, 29624, 58628,
                                                                       12968, 13103, 33728,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 63938, 0, 3,
                                                                       58628, 29840, 58988,
                                                                       13103, 13238, 33998,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 64388, 0, 3,
                                                                       58988, 30056, 59348,
                                                                       13238, 13373, 34268,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 64838, 0, 3,
                                                                       59348, 30272, 59708,
                                                                       13373, 13508, 34538,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 65288, 0, 3,
                                                                       60068, 31136, 60428,
                                                                       13778, 13913, 35348,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 65738, 0, 3,
                                                                       60428, 31352, 60788,
                                                                       13913, 14048, 35618,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 66188, 0, 3,
                                                                       60788, 31568, 61148,
                                                                       14048, 14183, 35888,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 66638, 0, 3,
                                                                       61148, 31784, 61508,
                                                                       14183, 14318, 36158,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 67088, 0, 3,
                                                                       61508, 32000, 61868,
                                                                       14318, 14453, 36428,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 67538, 0, 3,
                                                                       61868, 32216, 62228,
                                                                       14453, 14588, 36698,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 67988, 0, 3,
                                                                       62588, 33188, 63038,
                                                                       14858, 15023, 37628,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 68538, 0, 3,
                                                                       63038, 33458, 63488,
                                                                       15023, 15188, 37958,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 69088, 0, 3,
                                                                       63488, 33728, 63938,
                                                                       15188, 15353, 38288,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 69638, 0, 3,
                                                                       63938, 33998, 64388,
                                                                       15353, 15518, 38618,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 70188, 0, 3,
                                                                       64388, 34268, 64838,
                                                                       15518, 15683, 38948,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 70738, 0, 3,
                                                                       65288, 35348, 65738,
                                                                       16013, 16178, 39938,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 71288, 0, 3,
                                                                       65738, 35618, 66188,
                                                                       16178, 16343, 40268,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 71838, 0, 3,
                                                                       66188, 35888, 66638,
                                                                       16343, 16508, 40598,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 72388, 0, 3,
                                                                       66638, 36158, 67088,
                                                                       16508, 16673, 40928,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 72938, 0, 3,
                                                                       67088, 36428, 67538,
                                                                       16673, 16838, 41258,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73488, 3, 17168,
                                                                       17174, 41588, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73503, 3, 17174,
                                                                       17180, 41598, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73518, 3, 17180,
                                                                       17186, 41608, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73533, 3, 17186,
                                                                       17192, 41618, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73548, 3, 17192,
                                                                       17198, 41628, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73563, 3, 17198,
                                                                       17204, 41638, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73578, 3, 17204,
                                                                       17210, 41648, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73593, 3, 17210,
                                                                       17216, 41658, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73608, 3, 17216,
                                                                       17222, 41668, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73623, 3, 17222,
                                                                       17228, 41678, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73638, 3, 17228,
                                                                       17234, 41688, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73653, 3, 17234,
                                                                       17240, 41698, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73668, 3, 17240,
                                                                       17246, 41708, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73683, 3, 17246,
                                                                       17252, 41718, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73698, 3, 17264,
                                                                       17270, 41728, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73713, 3, 17270,
                                                                       17276, 41738, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73728, 3, 17276,
                                                                       17282, 41748, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73743, 3, 17282,
                                                                       17288, 41758, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73758, 3, 17288,
                                                                       17294, 41768, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73773, 3, 17294,
                                                                       17300, 41778, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73788, 3, 17300,
                                                                       17306, 41788, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73803, 3, 17306,
                                                                       17312, 41798, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73818, 3, 17312,
                                                                       17318, 41808, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73833, 3, 17318,
                                                                       17324, 41818, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73848, 3, 17324,
                                                                       17330, 41828, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73863, 3, 17330,
                                                                       17336, 41838, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73878, 3, 17336,
                                                                       17342, 41848, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 73893, 3, 17342,
                                                                       17348, 41858, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 73908, 0, 3,
                                                                       73488, 41588, 73503,
                                                                       17360, 17378, 41868,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 73953, 0, 3,
                                                                       73503, 41598, 73518,
                                                                       17378, 17396, 41898,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 73998, 0, 3,
                                                                       73518, 41608, 73533,
                                                                       17396, 17414, 41928,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74043, 0, 3,
                                                                       73533, 41618, 73548,
                                                                       17414, 17432, 41958,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74088, 0, 3,
                                                                       73548, 41628, 73563,
                                                                       17432, 17450, 41988,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74133, 0, 3,
                                                                       73563, 41638, 73578,
                                                                       17450, 17468, 42018,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74178, 0, 3,
                                                                       73578, 41648, 73593,
                                                                       17468, 17486, 42048,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74223, 0, 3,
                                                                       73593, 41658, 73608,
                                                                       17486, 17504, 42078,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74268, 0, 3,
                                                                       73608, 41668, 73623,
                                                                       17504, 17522, 42108,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74313, 0, 3,
                                                                       73623, 41678, 73638,
                                                                       17522, 17540, 42138,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74358, 0, 3,
                                                                       73638, 41688, 73653,
                                                                       17540, 17558, 42168,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74403, 0, 3,
                                                                       73653, 41698, 73668,
                                                                       17558, 17576, 42198,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74448, 0, 3,
                                                                       73668, 41708, 73683,
                                                                       17576, 17594, 42228,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74493, 0, 3,
                                                                       73698, 41728, 73713,
                                                                       17630, 17648, 42258,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74538, 0, 3,
                                                                       73713, 41738, 73728,
                                                                       17648, 17666, 42288,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74583, 0, 3,
                                                                       73728, 41748, 73743,
                                                                       17666, 17684, 42318,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74628, 0, 3,
                                                                       73743, 41758, 73758,
                                                                       17684, 17702, 42348,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74673, 0, 3,
                                                                       73758, 41768, 73773,
                                                                       17702, 17720, 42378,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74718, 0, 3,
                                                                       73773, 41778, 73788,
                                                                       17720, 17738, 42408,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74763, 0, 3,
                                                                       73788, 41788, 73803,
                                                                       17738, 17756, 42438,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74808, 0, 3,
                                                                       73803, 41798, 73818,
                                                                       17756, 17774, 42468,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74853, 0, 3,
                                                                       73818, 41808, 73833,
                                                                       17774, 17792, 42498,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74898, 0, 3,
                                                                       73833, 41818, 73848,
                                                                       17792, 17810, 42528,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74943, 0, 3,
                                                                       73848, 41828, 73863,
                                                                       17810, 17828, 42558,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 74988, 0, 3,
                                                                       73863, 41838, 73878,
                                                                       17828, 17846, 42588,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 75033, 0, 3,
                                                                       73878, 41848, 73893,
                                                                       17846, 17864, 42618,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 75078, 0, 3,
                                                                       73908, 41868, 73953,
                                                                       17900, 17936, 42648,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 75168, 0, 3,
                                                                       73953, 41898, 73998,
                                                                       17936, 17972, 42708,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 75258, 0, 3,
                                                                       73998, 41928, 74043,
                                                                       17972, 18008, 42768,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 75348, 0, 3,
                                                                       74043, 41958, 74088,
                                                                       18008, 18044, 42828,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 75438, 0, 3,
                                                                       74088, 41988, 74133,
                                                                       18044, 18080, 42888,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 75528, 0, 3,
                                                                       74133, 42018, 74178,
                                                                       18080, 18116, 42948,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 75618, 0, 3,
                                                                       74178, 42048, 74223,
                                                                       18116, 18152, 43008,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 75708, 0, 3,
                                                                       74223, 42078, 74268,
                                                                       18152, 18188, 43068,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 75798, 0, 3,
                                                                       74268, 42108, 74313,
                                                                       18188, 18224, 43128,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 75888, 0, 3,
                                                                       74313, 42138, 74358,
                                                                       18224, 18260, 43188,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 75978, 0, 3,
                                                                       74358, 42168, 74403,
                                                                       18260, 18296, 43248,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 76068, 0, 3,
                                                                       74403, 42198, 74448,
                                                                       18296, 18332, 43308,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 76158, 0, 3,
                                                                       74493, 42258, 74538,
                                                                       18404, 18440, 43368,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 76248, 0, 3,
                                                                       74538, 42288, 74583,
                                                                       18440, 18476, 43428,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 76338, 0, 3,
                                                                       74583, 42318, 74628,
                                                                       18476, 18512, 43488,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 76428, 0, 3,
                                                                       74628, 42348, 74673,
                                                                       18512, 18548, 43548,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 76518, 0, 3,
                                                                       74673, 42378, 74718,
                                                                       18548, 18584, 43608,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 76608, 0, 3,
                                                                       74718, 42408, 74763,
                                                                       18584, 18620, 43668,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 76698, 0, 3,
                                                                       74763, 42438, 74808,
                                                                       18620, 18656, 43728,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 76788, 0, 3,
                                                                       74808, 42468, 74853,
                                                                       18656, 18692, 43788,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 76878, 0, 3,
                                                                       74853, 42498, 74898,
                                                                       18692, 18728, 43848,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 76968, 0, 3,
                                                                       74898, 42528, 74943,
                                                                       18728, 18764, 43908,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77058, 0, 3,
                                                                       74943, 42558, 74988,
                                                                       18764, 18800, 43968,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77148, 0, 3,
                                                                       74988, 42588, 75033,
                                                                       18800, 18836, 44028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 77238, 0, 3,
                                                                       75078, 42648, 75168,
                                                                       18908, 18968, 44088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 77388, 0, 3,
                                                                       75168, 42708, 75258,
                                                                       18968, 19028, 44188,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 77538, 0, 3,
                                                                       75258, 42768, 75348,
                                                                       19028, 19088, 44288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 77688, 0, 3,
                                                                       75348, 42828, 75438,
                                                                       19088, 19148, 44388,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 77838, 0, 3,
                                                                       75438, 42888, 75528,
                                                                       19148, 19208, 44488,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 77988, 0, 3,
                                                                       75528, 42948, 75618,
                                                                       19208, 19268, 44588,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 78138, 0, 3,
                                                                       75618, 43008, 75708,
                                                                       19268, 19328, 44688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 78288, 0, 3,
                                                                       75708, 43068, 75798,
                                                                       19328, 19388, 44788,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 78438, 0, 3,
                                                                       75798, 43128, 75888,
                                                                       19388, 19448, 44888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 78588, 0, 3,
                                                                       75888, 43188, 75978,
                                                                       19448, 19508, 44988,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 78738, 0, 3,
                                                                       75978, 43248, 76068,
                                                                       19508, 19568, 45088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 78888, 0, 3,
                                                                       76158, 43368, 76248,
                                                                       19688, 19748, 45188,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79038, 0, 3,
                                                                       76248, 43428, 76338,
                                                                       19748, 19808, 45288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79188, 0, 3,
                                                                       76338, 43488, 76428,
                                                                       19808, 19868, 45388,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79338, 0, 3,
                                                                       76428, 43548, 76518,
                                                                       19868, 19928, 45488,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79488, 0, 3,
                                                                       76518, 43608, 76608,
                                                                       19928, 19988, 45588,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79638, 0, 3,
                                                                       76608, 43668, 76698,
                                                                       19988, 20048, 45688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79788, 0, 3,
                                                                       76698, 43728, 76788,
                                                                       20048, 20108, 45788,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79938, 0, 3,
                                                                       76788, 43788, 76878,
                                                                       20108, 20168, 45888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 80088, 0, 3,
                                                                       76878, 43848, 76968,
                                                                       20168, 20228, 45988,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 80238, 0, 3,
                                                                       76968, 43908, 77058,
                                                                       20228, 20288, 46088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 80388, 0, 3,
                                                                       77058, 43968, 77148,
                                                                       20288, 20348, 46188,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 80538, 0, 3,
                                                                       77238, 44088, 77388,
                                                                       20468, 20558, 46288,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 80763, 0, 3,
                                                                       77388, 44188, 77538,
                                                                       20558, 20648, 46438,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 80988, 0, 3,
                                                                       77538, 44288, 77688,
                                                                       20648, 20738, 46588,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 81213, 0, 3,
                                                                       77688, 44388, 77838,
                                                                       20738, 20828, 46738,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 81438, 0, 3,
                                                                       77838, 44488, 77988,
                                                                       20828, 20918, 46888,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 81663, 0, 3,
                                                                       77988, 44588, 78138,
                                                                       20918, 21008, 47038,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 81888, 0, 3,
                                                                       78138, 44688, 78288,
                                                                       21008, 21098, 47188,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 82113, 0, 3,
                                                                       78288, 44788, 78438,
                                                                       21098, 21188, 47338,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 82338, 0, 3,
                                                                       78438, 44888, 78588,
                                                                       21188, 21278, 47488,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 82563, 0, 3,
                                                                       78588, 44988, 78738,
                                                                       21278, 21368, 47638,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 82788, 0, 3,
                                                                       78888, 45188, 79038,
                                                                       21548, 21638, 47788,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 83013, 0, 3,
                                                                       79038, 45288, 79188,
                                                                       21638, 21728, 47938,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 83238, 0, 3,
                                                                       79188, 45388, 79338,
                                                                       21728, 21818, 48088,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 83463, 0, 3,
                                                                       79338, 45488, 79488,
                                                                       21818, 21908, 48238,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 83688, 0, 3,
                                                                       79488, 45588, 79638,
                                                                       21908, 21998, 48388,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 83913, 0, 3,
                                                                       79638, 45688, 79788,
                                                                       21998, 22088, 48538,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 84138, 0, 3,
                                                                       79788, 45788, 79938,
                                                                       22088, 22178, 48688,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 84363, 0, 3,
                                                                       79938, 45888, 80088,
                                                                       22178, 22268, 48838,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 84588, 0, 3,
                                                                       80088, 45988, 80238,
                                                                       22268, 22358, 48988,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 84813, 0, 3,
                                                                       80238, 46088, 80388,
                                                                       22358, 22448, 49138,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 85038, 0, 3,
                                                                       80538, 46288, 80763,
                                                                       22628, 22754, 49288,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 85353, 0, 3,
                                                                       80763, 46438, 80988,
                                                                       22754, 22880, 49498,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 85668, 0, 3,
                                                                       80988, 46588, 81213,
                                                                       22880, 23006, 49708,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 85983, 0, 3,
                                                                       81213, 46738, 81438,
                                                                       23006, 23132, 49918,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 86298, 0, 3,
                                                                       81438, 46888, 81663,
                                                                       23132, 23258, 50128,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 86613, 0, 3,
                                                                       81663, 47038, 81888,
                                                                       23258, 23384, 50338,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 86928, 0, 3,
                                                                       81888, 47188, 82113,
                                                                       23384, 23510, 50548,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 87243, 0, 3,
                                                                       82113, 47338, 82338,
                                                                       23510, 23636, 50758,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 87558, 0, 3,
                                                                       82338, 47488, 82563,
                                                                       23636, 23762, 50968,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 87873, 0, 3,
                                                                       82788, 47788, 83013,
                                                                       24014, 24140, 51178,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 88188, 0, 3,
                                                                       83013, 47938, 83238,
                                                                       24140, 24266, 51388,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 88503, 0, 3,
                                                                       83238, 48088, 83463,
                                                                       24266, 24392, 51598,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 88818, 0, 3,
                                                                       83463, 48238, 83688,
                                                                       24392, 24518, 51808,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 89133, 0, 3,
                                                                       83688, 48388, 83913,
                                                                       24518, 24644, 52018,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 89448, 0, 3,
                                                                       83913, 48538, 84138,
                                                                       24644, 24770, 52228,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 89763, 0, 3,
                                                                       84138, 48688, 84363,
                                                                       24770, 24896, 52438,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 90078, 0, 3,
                                                                       84363, 48838, 84588,
                                                                       24896, 25022, 52648,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 90393, 0, 3,
                                                                       84588, 48988, 84813,
                                                                       25022, 25148, 52858,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 90708, 0, 3,
                                                                       85038, 49288, 85353,
                                                                       25400, 25568, 53068,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 91128, 0, 3,
                                                                       85353, 49498, 85668,
                                                                       25568, 25736, 53348,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 91548, 0, 3,
                                                                       85668, 49708, 85983,
                                                                       25736, 25904, 53628,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 91968, 0, 3,
                                                                       85983, 49918, 86298,
                                                                       25904, 26072, 53908,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 92388, 0, 3,
                                                                       86298, 50128, 86613,
                                                                       26072, 26240, 54188,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 92808, 0, 3,
                                                                       86613, 50338, 86928,
                                                                       26240, 26408, 54468,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 93228, 0, 3,
                                                                       86928, 50548, 87243,
                                                                       26408, 26576, 54748,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 93648, 0, 3,
                                                                       87243, 50758, 87558,
                                                                       26576, 26744, 55028,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 94068, 0, 3,
                                                                       87873, 51178, 88188,
                                                                       27080, 27248, 55308,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 94488, 0, 3,
                                                                       88188, 51388, 88503,
                                                                       27248, 27416, 55588,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 94908, 0, 3,
                                                                       88503, 51598, 88818,
                                                                       27416, 27584, 55868,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 95328, 0, 3,
                                                                       88818, 51808, 89133,
                                                                       27584, 27752, 56148,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 95748, 0, 3,
                                                                       89133, 52018, 89448,
                                                                       27752, 27920, 56428,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 96168, 0, 3,
                                                                       89448, 52228, 89763,
                                                                       27920, 28088, 56708,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 96588, 0, 3,
                                                                       89763, 52438, 90078,
                                                                       28088, 28256, 56988,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 97008, 0, 3,
                                                                       90078, 52648, 90393,
                                                                       28256, 28424, 57268,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 97428, 0, 3,
                                                                       90708, 53068, 91128,
                                                                       28760, 28976, 57548,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 97968, 0, 3,
                                                                       91128, 53348, 91548,
                                                                       28976, 29192, 57908,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 98508, 0, 3,
                                                                       91548, 53628, 91968,
                                                                       29192, 29408, 58268,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 99048, 0, 3,
                                                                       91968, 53908, 92388,
                                                                       29408, 29624, 58628,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 99588, 0, 3,
                                                                       92388, 54188, 92808,
                                                                       29624, 29840, 58988,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 100128, 0, 3,
                                                                       92808, 54468, 93228,
                                                                       29840, 30056, 59348,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 100668, 0, 3,
                                                                       93228, 54748, 93648,
                                                                       30056, 30272, 59708,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 101208, 0, 3,
                                                                       94068, 55308, 94488,
                                                                       30704, 30920, 60068,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 101748, 0, 3,
                                                                       94488, 55588, 94908,
                                                                       30920, 31136, 60428,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 102288, 0, 3,
                                                                       94908, 55868, 95328,
                                                                       31136, 31352, 60788,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 102828, 0, 3,
                                                                       95328, 56148, 95748,
                                                                       31352, 31568, 61148,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 103368, 0, 3,
                                                                       95748, 56428, 96168,
                                                                       31568, 31784, 61508,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 103908, 0, 3,
                                                                       96168, 56708, 96588,
                                                                       31784, 32000, 61868,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 104448, 0, 3,
                                                                       96588, 56988, 97008,
                                                                       32000, 32216, 62228,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 104988, 0, 3,
                                                                       97428, 57548, 97968,
                                                                       32648, 32918, 62588,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 105663, 0, 3,
                                                                       97968, 57908, 98508,
                                                                       32918, 33188, 63038,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 106338, 0, 3,
                                                                       98508, 58268, 99048,
                                                                       33188, 33458, 63488,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 107013, 0, 3,
                                                                       99048, 58628, 99588,
                                                                       33458, 33728, 63938,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 107688, 0, 3,
                                                                       99588, 58988, 100128,
                                                                       33728, 33998, 64388,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 108363, 0, 3,
                                                                       100128, 59348, 100668,
                                                                       33998, 34268, 64838,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 109038, 0, 3,
                                                                       101208, 60068, 101748,
                                                                       34808, 35078, 65288,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 109713, 0, 3,
                                                                       101748, 60428, 102288,
                                                                       35078, 35348, 65738,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 110388, 0, 3,
                                                                       102288, 60788, 102828,
                                                                       35348, 35618, 66188,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 111063, 0, 3,
                                                                       102828, 61148, 103368,
                                                                       35618, 35888, 66638,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 111738, 0, 3,
                                                                       103368, 61508, 103908,
                                                                       35888, 36158, 67088,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 112413, 0, 3,
                                                                       103908, 61868, 104448,
                                                                       36158, 36428, 67538,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 113088, 0, 3,
                                                                       104988, 62588, 105663,
                                                                       36968, 37298, 67988,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 113913, 0, 3,
                                                                       105663, 63038, 106338,
                                                                       37298, 37628, 68538,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 114738, 0, 3,
                                                                       106338, 63488, 107013,
                                                                       37628, 37958, 69088,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 115563, 0, 3,
                                                                       107013, 63938, 107688,
                                                                       37958, 38288, 69638,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 116388, 0, 3,
                                                                       107688, 64388, 108363,
                                                                       38288, 38618, 70188,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 117213, 0, 3,
                                                                       109038, 65288, 109713,
                                                                       39278, 39608, 70738,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 118038, 0, 3,
                                                                       109713, 65738, 110388,
                                                                       39608, 39938, 71288,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 118863, 0, 3,
                                                                       110388, 66188, 111063,
                                                                       39938, 40268, 71838,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 119688, 0, 3,
                                                                       111063, 66638, 111738,
                                                                       40268, 40598, 72388,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 120513, 0, 3,
                                                                       111738, 67088, 112413,
                                                                       40598, 40928, 72938,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121338, 3, 41588,
                                                                       41598, 73518, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121359, 3, 41598,
                                                                       41608, 73533, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121380, 3, 41608,
                                                                       41618, 73548, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121401, 3, 41618,
                                                                       41628, 73563, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121422, 3, 41628,
                                                                       41638, 73578, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121443, 3, 41638,
                                                                       41648, 73593, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121464, 3, 41648,
                                                                       41658, 73608, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121485, 3, 41658,
                                                                       41668, 73623, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121506, 3, 41668,
                                                                       41678, 73638, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121527, 3, 41678,
                                                                       41688, 73653, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121548, 3, 41688,
                                                                       41698, 73668, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121569, 3, 41698,
                                                                       41708, 73683, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121590, 3, 41728,
                                                                       41738, 73728, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121611, 3, 41738,
                                                                       41748, 73743, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121632, 3, 41748,
                                                                       41758, 73758, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121653, 3, 41758,
                                                                       41768, 73773, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121674, 3, 41768,
                                                                       41778, 73788, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121695, 3, 41778,
                                                                       41788, 73803, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121716, 3, 41788,
                                                                       41798, 73818, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121737, 3, 41798,
                                                                       41808, 73833, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121758, 3, 41808,
                                                                       41818, 73848, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121779, 3, 41818,
                                                                       41828, 73863, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121800, 3, 41828,
                                                                       41838, 73878, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 121821, 3, 41838,
                                                                       41848, 73893, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 121842, 0, 3,
                                                                       121338, 73518, 121359,
                                                                       41868, 41898, 73998,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 121905, 0, 3,
                                                                       121359, 73533, 121380,
                                                                       41898, 41928, 74043,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 121968, 0, 3,
                                                                       121380, 73548, 121401,
                                                                       41928, 41958, 74088,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 122031, 0, 3,
                                                                       121401, 73563, 121422,
                                                                       41958, 41988, 74133,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 122094, 0, 3,
                                                                       121422, 73578, 121443,
                                                                       41988, 42018, 74178,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 122157, 0, 3,
                                                                       121443, 73593, 121464,
                                                                       42018, 42048, 74223,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 122220, 0, 3,
                                                                       121464, 73608, 121485,
                                                                       42048, 42078, 74268,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 122283, 0, 3,
                                                                       121485, 73623, 121506,
                                                                       42078, 42108, 74313,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 122346, 0, 3,
                                                                       121506, 73638, 121527,
                                                                       42108, 42138, 74358,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 122409, 0, 3,
                                                                       121527, 73653, 121548,
                                                                       42138, 42168, 74403,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 122472, 0, 3,
                                                                       121548, 73668, 121569,
                                                                       42168, 42198, 74448,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 122535, 0, 3,
                                                                       121590, 73728, 121611,
                                                                       42258, 42288, 74583,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 122598, 0, 3,
                                                                       121611, 73743, 121632,
                                                                       42288, 42318, 74628,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 122661, 0, 3,
                                                                       121632, 73758, 121653,
                                                                       42318, 42348, 74673,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 122724, 0, 3,
                                                                       121653, 73773, 121674,
                                                                       42348, 42378, 74718,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 122787, 0, 3,
                                                                       121674, 73788, 121695,
                                                                       42378, 42408, 74763,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 122850, 0, 3,
                                                                       121695, 73803, 121716,
                                                                       42408, 42438, 74808,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 122913, 0, 3,
                                                                       121716, 73818, 121737,
                                                                       42438, 42468, 74853,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 122976, 0, 3,
                                                                       121737, 73833, 121758,
                                                                       42468, 42498, 74898,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 123039, 0, 3,
                                                                       121758, 73848, 121779,
                                                                       42498, 42528, 74943,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 123102, 0, 3,
                                                                       121779, 73863, 121800,
                                                                       42528, 42558, 74988,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 123165, 0, 3,
                                                                       121800, 73878, 121821,
                                                                       42558, 42588, 75033,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 123228, 0, 3,
                                                                       121842, 73998, 121905,
                                                                       42648, 42708, 75258,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 123354, 0, 3,
                                                                       121905, 74043, 121968,
                                                                       42708, 42768, 75348,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 123480, 0, 3,
                                                                       121968, 74088, 122031,
                                                                       42768, 42828, 75438,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 123606, 0, 3,
                                                                       122031, 74133, 122094,
                                                                       42828, 42888, 75528,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 123732, 0, 3,
                                                                       122094, 74178, 122157,
                                                                       42888, 42948, 75618,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 123858, 0, 3,
                                                                       122157, 74223, 122220,
                                                                       42948, 43008, 75708,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 123984, 0, 3,
                                                                       122220, 74268, 122283,
                                                                       43008, 43068, 75798,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 124110, 0, 3,
                                                                       122283, 74313, 122346,
                                                                       43068, 43128, 75888,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 124236, 0, 3,
                                                                       122346, 74358, 122409,
                                                                       43128, 43188, 75978,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 124362, 0, 3,
                                                                       122409, 74403, 122472,
                                                                       43188, 43248, 76068,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 124488, 0, 3,
                                                                       122535, 74583, 122598,
                                                                       43368, 43428, 76338,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 124614, 0, 3,
                                                                       122598, 74628, 122661,
                                                                       43428, 43488, 76428,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 124740, 0, 3,
                                                                       122661, 74673, 122724,
                                                                       43488, 43548, 76518,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 124866, 0, 3,
                                                                       122724, 74718, 122787,
                                                                       43548, 43608, 76608,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 124992, 0, 3,
                                                                       122787, 74763, 122850,
                                                                       43608, 43668, 76698,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 125118, 0, 3,
                                                                       122850, 74808, 122913,
                                                                       43668, 43728, 76788,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 125244, 0, 3,
                                                                       122913, 74853, 122976,
                                                                       43728, 43788, 76878,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 125370, 0, 3,
                                                                       122976, 74898, 123039,
                                                                       43788, 43848, 76968,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 125496, 0, 3,
                                                                       123039, 74943, 123102,
                                                                       43848, 43908, 77058,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 125622, 0, 3,
                                                                       123102, 74988, 123165,
                                                                       43908, 43968, 77148,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 125748, 0, 3,
                                                                       123228, 75258, 123354,
                                                                       44088, 44188, 77538,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 125958, 0, 3,
                                                                       123354, 75348, 123480,
                                                                       44188, 44288, 77688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 126168, 0, 3,
                                                                       123480, 75438, 123606,
                                                                       44288, 44388, 77838,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 126378, 0, 3,
                                                                       123606, 75528, 123732,
                                                                       44388, 44488, 77988,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 126588, 0, 3,
                                                                       123732, 75618, 123858,
                                                                       44488, 44588, 78138,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 126798, 0, 3,
                                                                       123858, 75708, 123984,
                                                                       44588, 44688, 78288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 127008, 0, 3,
                                                                       123984, 75798, 124110,
                                                                       44688, 44788, 78438,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 127218, 0, 3,
                                                                       124110, 75888, 124236,
                                                                       44788, 44888, 78588,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 127428, 0, 3,
                                                                       124236, 75978, 124362,
                                                                       44888, 44988, 78738,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 127638, 0, 3,
                                                                       124488, 76338, 124614,
                                                                       45188, 45288, 79188,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 127848, 0, 3,
                                                                       124614, 76428, 124740,
                                                                       45288, 45388, 79338,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 128058, 0, 3,
                                                                       124740, 76518, 124866,
                                                                       45388, 45488, 79488,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 128268, 0, 3,
                                                                       124866, 76608, 124992,
                                                                       45488, 45588, 79638,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 128478, 0, 3,
                                                                       124992, 76698, 125118,
                                                                       45588, 45688, 79788,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 128688, 0, 3,
                                                                       125118, 76788, 125244,
                                                                       45688, 45788, 79938,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 128898, 0, 3,
                                                                       125244, 76878, 125370,
                                                                       45788, 45888, 80088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 129108, 0, 3,
                                                                       125370, 76968, 125496,
                                                                       45888, 45988, 80238,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 129318, 0, 3,
                                                                       125496, 77058, 125622,
                                                                       45988, 46088, 80388,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 129528, 0, 3,
                                                                       125748, 77538, 125958,
                                                                       46288, 46438, 80988,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 129843, 0, 3,
                                                                       125958, 77688, 126168,
                                                                       46438, 46588, 81213,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 130158, 0, 3,
                                                                       126168, 77838, 126378,
                                                                       46588, 46738, 81438,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 130473, 0, 3,
                                                                       126378, 77988, 126588,
                                                                       46738, 46888, 81663,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 130788, 0, 3,
                                                                       126588, 78138, 126798,
                                                                       46888, 47038, 81888,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 131103, 0, 3,
                                                                       126798, 78288, 127008,
                                                                       47038, 47188, 82113,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 131418, 0, 3,
                                                                       127008, 78438, 127218,
                                                                       47188, 47338, 82338,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 131733, 0, 3,
                                                                       127218, 78588, 127428,
                                                                       47338, 47488, 82563,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 132048, 0, 3,
                                                                       127638, 79188, 127848,
                                                                       47788, 47938, 83238,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 132363, 0, 3,
                                                                       127848, 79338, 128058,
                                                                       47938, 48088, 83463,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 132678, 0, 3,
                                                                       128058, 79488, 128268,
                                                                       48088, 48238, 83688,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 132993, 0, 3,
                                                                       128268, 79638, 128478,
                                                                       48238, 48388, 83913,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 133308, 0, 3,
                                                                       128478, 79788, 128688,
                                                                       48388, 48538, 84138,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 133623, 0, 3,
                                                                       128688, 79938, 128898,
                                                                       48538, 48688, 84363,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 133938, 0, 3,
                                                                       128898, 80088, 129108,
                                                                       48688, 48838, 84588,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 134253, 0, 3,
                                                                       129108, 80238, 129318,
                                                                       48838, 48988, 84813,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 134568, 0, 3,
                                                                       129528, 80988, 129843,
                                                                       49288, 49498, 85668,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 135009, 0, 3,
                                                                       129843, 81213, 130158,
                                                                       49498, 49708, 85983,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 135450, 0, 3,
                                                                       130158, 81438, 130473,
                                                                       49708, 49918, 86298,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 135891, 0, 3,
                                                                       130473, 81663, 130788,
                                                                       49918, 50128, 86613,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 136332, 0, 3,
                                                                       130788, 81888, 131103,
                                                                       50128, 50338, 86928,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 136773, 0, 3,
                                                                       131103, 82113, 131418,
                                                                       50338, 50548, 87243,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 137214, 0, 3,
                                                                       131418, 82338, 131733,
                                                                       50548, 50758, 87558,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 137655, 0, 3,
                                                                       132048, 83238, 132363,
                                                                       51178, 51388, 88503,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 138096, 0, 3,
                                                                       132363, 83463, 132678,
                                                                       51388, 51598, 88818,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 138537, 0, 3,
                                                                       132678, 83688, 132993,
                                                                       51598, 51808, 89133,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 138978, 0, 3,
                                                                       132993, 83913, 133308,
                                                                       51808, 52018, 89448,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 139419, 0, 3,
                                                                       133308, 84138, 133623,
                                                                       52018, 52228, 89763,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 139860, 0, 3,
                                                                       133623, 84363, 133938,
                                                                       52228, 52438, 90078,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 140301, 0, 3,
                                                                       133938, 84588, 134253,
                                                                       52438, 52648, 90393,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 140742, 0, 3,
                                                                       134568, 85668, 135009,
                                                                       53068, 53348, 91548,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 141330, 0, 3,
                                                                       135009, 85983, 135450,
                                                                       53348, 53628, 91968,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 141918, 0, 3,
                                                                       135450, 86298, 135891,
                                                                       53628, 53908, 92388,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 142506, 0, 3,
                                                                       135891, 86613, 136332,
                                                                       53908, 54188, 92808,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 143094, 0, 3,
                                                                       136332, 86928, 136773,
                                                                       54188, 54468, 93228,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 143682, 0, 3,
                                                                       136773, 87243, 137214,
                                                                       54468, 54748, 93648,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 144270, 0, 3,
                                                                       137655, 88503, 138096,
                                                                       55308, 55588, 94908,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 144858, 0, 3,
                                                                       138096, 88818, 138537,
                                                                       55588, 55868, 95328,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 145446, 0, 3,
                                                                       138537, 89133, 138978,
                                                                       55868, 56148, 95748,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 146034, 0, 3,
                                                                       138978, 89448, 139419,
                                                                       56148, 56428, 96168,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 146622, 0, 3,
                                                                       139419, 89763, 139860,
                                                                       56428, 56708, 96588,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 147210, 0, 3,
                                                                       139860, 90078, 140301,
                                                                       56708, 56988, 97008,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 147798, 0, 3,
                                                                       140742, 91548, 141330,
                                                                       57548, 57908, 98508,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 148554, 0, 3,
                                                                       141330, 91968, 141918,
                                                                       57908, 58268, 99048,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 149310, 0, 3,
                                                                       141918, 92388, 142506,
                                                                       58268, 58628, 99588,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 150066, 0, 3,
                                                                       142506, 92808, 143094,
                                                                       58628, 58988, 100128,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 150822, 0, 3,
                                                                       143094, 93228, 143682,
                                                                       58988, 59348, 100668,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 151578, 0, 3,
                                                                       144270, 94908, 144858,
                                                                       60068, 60428, 102288,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 152334, 0, 3,
                                                                       144858, 95328, 145446,
                                                                       60428, 60788, 102828,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 153090, 0, 3,
                                                                       145446, 95748, 146034,
                                                                       60788, 61148, 103368,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 153846, 0, 3,
                                                                       146034, 96168, 146622,
                                                                       61148, 61508, 103908,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 154602, 0, 3,
                                                                       146622, 96588, 147210,
                                                                       61508, 61868, 104448,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 155358, 0, 3,
                                                                       147798, 98508, 148554,
                                                                       62588, 63038, 106338,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 156303, 0, 3,
                                                                       148554, 99048, 149310,
                                                                       63038, 63488, 107013,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 157248, 0, 3,
                                                                       149310, 99588, 150066,
                                                                       63488, 63938, 107688,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 158193, 0, 3,
                                                                       150066, 100128, 150822,
                                                                       63938, 64388, 108363,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 159138, 0, 3,
                                                                       151578, 102288, 152334,
                                                                       65288, 65738, 110388,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 160083, 0, 3,
                                                                       152334, 102828, 153090,
                                                                       65738, 66188, 111063,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 161028, 0, 3,
                                                                       153090, 103368, 153846,
                                                                       66188, 66638, 111738,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 161973, 0, 3,
                                                                       153846, 103908, 154602,
                                                                       66638, 67088, 112413,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 162918, 0, 3,
                                                                       155358, 106338, 156303,
                                                                       67988, 68538, 114738,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 164073, 0, 3,
                                                                       156303, 107013, 157248,
                                                                       68538, 69088, 115563,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 165228, 0, 3,
                                                                       157248, 107688, 158193,
                                                                       69088, 69638, 116388,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 166383, 0, 3,
                                                                       159138, 110388, 160083,
                                                                       70738, 71288, 118863,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 167538, 0, 3,
                                                                       160083, 111063, 161028,
                                                                       71288, 71838, 119688,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 168693, 0, 3,
                                                                       161028, 111738, 161973,
                                                                       71838, 72388, 120513,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 169848, 3, 73488,
                                                                       73503, 121338, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 169876, 3, 73503,
                                                                       73518, 121359, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 169904, 3, 73518,
                                                                       73533, 121380, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 169932, 3, 73533,
                                                                       73548, 121401, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 169960, 3, 73548,
                                                                       73563, 121422, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 169988, 3, 73563,
                                                                       73578, 121443, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 170016, 3, 73578,
                                                                       73593, 121464, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 170044, 3, 73593,
                                                                       73608, 121485, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 170072, 3, 73608,
                                                                       73623, 121506, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 170100, 3, 73623,
                                                                       73638, 121527, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 170128, 3, 73638,
                                                                       73653, 121548, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 170156, 3, 73653,
                                                                       73668, 121569, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 170184, 3, 73698,
                                                                       73713, 121590, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 170212, 3, 73713,
                                                                       73728, 121611, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 170240, 3, 73728,
                                                                       73743, 121632, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 170268, 3, 73743,
                                                                       73758, 121653, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 170296, 3, 73758,
                                                                       73773, 121674, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 170324, 3, 73773,
                                                                       73788, 121695, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 170352, 3, 73788,
                                                                       73803, 121716, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 170380, 3, 73803,
                                                                       73818, 121737, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 170408, 3, 73818,
                                                                       73833, 121758, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 170436, 3, 73833,
                                                                       73848, 121779, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 170464, 3, 73848,
                                                                       73863, 121800, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 170492, 3, 73863,
                                                                       73878, 121821, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 170520, 0, 3,
                                                                       169848, 121338, 169876,
                                                                       73908, 73953, 121842,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 170604, 0, 3,
                                                                       169876, 121359, 169904,
                                                                       73953, 73998, 121905,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 170688, 0, 3,
                                                                       169904, 121380, 169932,
                                                                       73998, 74043, 121968,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 170772, 0, 3,
                                                                       169932, 121401, 169960,
                                                                       74043, 74088, 122031,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 170856, 0, 3,
                                                                       169960, 121422, 169988,
                                                                       74088, 74133, 122094,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 170940, 0, 3,
                                                                       169988, 121443, 170016,
                                                                       74133, 74178, 122157,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 171024, 0, 3,
                                                                       170016, 121464, 170044,
                                                                       74178, 74223, 122220,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 171108, 0, 3,
                                                                       170044, 121485, 170072,
                                                                       74223, 74268, 122283,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 171192, 0, 3,
                                                                       170072, 121506, 170100,
                                                                       74268, 74313, 122346,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 171276, 0, 3,
                                                                       170100, 121527, 170128,
                                                                       74313, 74358, 122409,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 171360, 0, 3,
                                                                       170128, 121548, 170156,
                                                                       74358, 74403, 122472,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 171444, 0, 3,
                                                                       170184, 121590, 170212,
                                                                       74493, 74538, 122535,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 171528, 0, 3,
                                                                       170212, 121611, 170240,
                                                                       74538, 74583, 122598,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 171612, 0, 3,
                                                                       170240, 121632, 170268,
                                                                       74583, 74628, 122661,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 171696, 0, 3,
                                                                       170268, 121653, 170296,
                                                                       74628, 74673, 122724,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 171780, 0, 3,
                                                                       170296, 121674, 170324,
                                                                       74673, 74718, 122787,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 171864, 0, 3,
                                                                       170324, 121695, 170352,
                                                                       74718, 74763, 122850,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 171948, 0, 3,
                                                                       170352, 121716, 170380,
                                                                       74763, 74808, 122913,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 172032, 0, 3,
                                                                       170380, 121737, 170408,
                                                                       74808, 74853, 122976,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 172116, 0, 3,
                                                                       170408, 121758, 170436,
                                                                       74853, 74898, 123039,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 172200, 0, 3,
                                                                       170436, 121779, 170464,
                                                                       74898, 74943, 123102,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 172284, 0, 3,
                                                                       170464, 121800, 170492,
                                                                       74943, 74988, 123165,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 172368, 0, 3,
                                                                       170520, 121842, 170604,
                                                                       75078, 75168, 123228,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 172536, 0, 3,
                                                                       170604, 121905, 170688,
                                                                       75168, 75258, 123354,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 172704, 0, 3,
                                                                       170688, 121968, 170772,
                                                                       75258, 75348, 123480,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 172872, 0, 3,
                                                                       170772, 122031, 170856,
                                                                       75348, 75438, 123606,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 173040, 0, 3,
                                                                       170856, 122094, 170940,
                                                                       75438, 75528, 123732,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 173208, 0, 3,
                                                                       170940, 122157, 171024,
                                                                       75528, 75618, 123858,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 173376, 0, 3,
                                                                       171024, 122220, 171108,
                                                                       75618, 75708, 123984,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 173544, 0, 3,
                                                                       171108, 122283, 171192,
                                                                       75708, 75798, 124110,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 173712, 0, 3,
                                                                       171192, 122346, 171276,
                                                                       75798, 75888, 124236,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 173880, 0, 3,
                                                                       171276, 122409, 171360,
                                                                       75888, 75978, 124362,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 174048, 0, 3,
                                                                       171444, 122535, 171528,
                                                                       76158, 76248, 124488,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 174216, 0, 3,
                                                                       171528, 122598, 171612,
                                                                       76248, 76338, 124614,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 174384, 0, 3,
                                                                       171612, 122661, 171696,
                                                                       76338, 76428, 124740,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 174552, 0, 3,
                                                                       171696, 122724, 171780,
                                                                       76428, 76518, 124866,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 174720, 0, 3,
                                                                       171780, 122787, 171864,
                                                                       76518, 76608, 124992,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 174888, 0, 3,
                                                                       171864, 122850, 171948,
                                                                       76608, 76698, 125118,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 175056, 0, 3,
                                                                       171948, 122913, 172032,
                                                                       76698, 76788, 125244,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 175224, 0, 3,
                                                                       172032, 122976, 172116,
                                                                       76788, 76878, 125370,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 175392, 0, 3,
                                                                       172116, 123039, 172200,
                                                                       76878, 76968, 125496,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 175560, 0, 3,
                                                                       172200, 123102, 172284,
                                                                       76968, 77058, 125622,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 175728, 0, 3,
                                                                       172368, 123228, 172536,
                                                                       77238, 77388, 125748,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 176008, 0, 3,
                                                                       172536, 123354, 172704,
                                                                       77388, 77538, 125958,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 176288, 0, 3,
                                                                       172704, 123480, 172872,
                                                                       77538, 77688, 126168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 176568, 0, 3,
                                                                       172872, 123606, 173040,
                                                                       77688, 77838, 126378,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 176848, 0, 3,
                                                                       173040, 123732, 173208,
                                                                       77838, 77988, 126588,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 177128, 0, 3,
                                                                       173208, 123858, 173376,
                                                                       77988, 78138, 126798,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 177408, 0, 3,
                                                                       173376, 123984, 173544,
                                                                       78138, 78288, 127008,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 177688, 0, 3,
                                                                       173544, 124110, 173712,
                                                                       78288, 78438, 127218,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 177968, 0, 3,
                                                                       173712, 124236, 173880,
                                                                       78438, 78588, 127428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 178248, 0, 3,
                                                                       174048, 124488, 174216,
                                                                       78888, 79038, 127638,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 178528, 0, 3,
                                                                       174216, 124614, 174384,
                                                                       79038, 79188, 127848,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 178808, 0, 3,
                                                                       174384, 124740, 174552,
                                                                       79188, 79338, 128058,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 179088, 0, 3,
                                                                       174552, 124866, 174720,
                                                                       79338, 79488, 128268,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 179368, 0, 3,
                                                                       174720, 124992, 174888,
                                                                       79488, 79638, 128478,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 179648, 0, 3,
                                                                       174888, 125118, 175056,
                                                                       79638, 79788, 128688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 179928, 0, 3,
                                                                       175056, 125244, 175224,
                                                                       79788, 79938, 128898,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 180208, 0, 3,
                                                                       175224, 125370, 175392,
                                                                       79938, 80088, 129108,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 180488, 0, 3,
                                                                       175392, 125496, 175560,
                                                                       80088, 80238, 129318,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 180768, 0, 3,
                                                                       175728, 125748, 176008,
                                                                       80538, 80763, 129528,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 181188, 0, 3,
                                                                       176008, 125958, 176288,
                                                                       80763, 80988, 129843,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 181608, 0, 3,
                                                                       176288, 126168, 176568,
                                                                       80988, 81213, 130158,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 182028, 0, 3,
                                                                       176568, 126378, 176848,
                                                                       81213, 81438, 130473,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 182448, 0, 3,
                                                                       176848, 126588, 177128,
                                                                       81438, 81663, 130788,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 182868, 0, 3,
                                                                       177128, 126798, 177408,
                                                                       81663, 81888, 131103,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 183288, 0, 3,
                                                                       177408, 127008, 177688,
                                                                       81888, 82113, 131418,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 183708, 0, 3,
                                                                       177688, 127218, 177968,
                                                                       82113, 82338, 131733,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 184128, 0, 3,
                                                                       178248, 127638, 178528,
                                                                       82788, 83013, 132048,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 184548, 0, 3,
                                                                       178528, 127848, 178808,
                                                                       83013, 83238, 132363,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 184968, 0, 3,
                                                                       178808, 128058, 179088,
                                                                       83238, 83463, 132678,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 185388, 0, 3,
                                                                       179088, 128268, 179368,
                                                                       83463, 83688, 132993,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 185808, 0, 3,
                                                                       179368, 128478, 179648,
                                                                       83688, 83913, 133308,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 186228, 0, 3,
                                                                       179648, 128688, 179928,
                                                                       83913, 84138, 133623,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 186648, 0, 3,
                                                                       179928, 128898, 180208,
                                                                       84138, 84363, 133938,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 187068, 0, 3,
                                                                       180208, 129108, 180488,
                                                                       84363, 84588, 134253,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 187488, 0, 3,
                                                                       180768, 129528, 181188,
                                                                       85038, 85353, 134568,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 188076, 0, 3,
                                                                       181188, 129843, 181608,
                                                                       85353, 85668, 135009,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 188664, 0, 3,
                                                                       181608, 130158, 182028,
                                                                       85668, 85983, 135450,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 189252, 0, 3,
                                                                       182028, 130473, 182448,
                                                                       85983, 86298, 135891,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 189840, 0, 3,
                                                                       182448, 130788, 182868,
                                                                       86298, 86613, 136332,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 190428, 0, 3,
                                                                       182868, 131103, 183288,
                                                                       86613, 86928, 136773,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 191016, 0, 3,
                                                                       183288, 131418, 183708,
                                                                       86928, 87243, 137214,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 191604, 0, 3,
                                                                       184128, 132048, 184548,
                                                                       87873, 88188, 137655,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 192192, 0, 3,
                                                                       184548, 132363, 184968,
                                                                       88188, 88503, 138096,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 192780, 0, 3,
                                                                       184968, 132678, 185388,
                                                                       88503, 88818, 138537,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 193368, 0, 3,
                                                                       185388, 132993, 185808,
                                                                       88818, 89133, 138978,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 193956, 0, 3,
                                                                       185808, 133308, 186228,
                                                                       89133, 89448, 139419,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 194544, 0, 3,
                                                                       186228, 133623, 186648,
                                                                       89448, 89763, 139860,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 195132, 0, 3,
                                                                       186648, 133938, 187068,
                                                                       89763, 90078, 140301,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 195720, 0, 3,
                                                                       187488, 134568, 188076,
                                                                       90708, 91128, 140742,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 196504, 0, 3,
                                                                       188076, 135009, 188664,
                                                                       91128, 91548, 141330,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 197288, 0, 3,
                                                                       188664, 135450, 189252,
                                                                       91548, 91968, 141918,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 198072, 0, 3,
                                                                       189252, 135891, 189840,
                                                                       91968, 92388, 142506,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 198856, 0, 3,
                                                                       189840, 136332, 190428,
                                                                       92388, 92808, 143094,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 199640, 0, 3,
                                                                       190428, 136773, 191016,
                                                                       92808, 93228, 143682,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 200424, 0, 3,
                                                                       191604, 137655, 192192,
                                                                       94068, 94488, 144270,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 201208, 0, 3,
                                                                       192192, 138096, 192780,
                                                                       94488, 94908, 144858,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 201992, 0, 3,
                                                                       192780, 138537, 193368,
                                                                       94908, 95328, 145446,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 202776, 0, 3,
                                                                       193368, 138978, 193956,
                                                                       95328, 95748, 146034,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 203560, 0, 3,
                                                                       193956, 139419, 194544,
                                                                       95748, 96168, 146622,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 204344, 0, 3,
                                                                       194544, 139860, 195132,
                                                                       96168, 96588, 147210,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 205128, 0, 3,
                                                                       195720, 140742, 196504,
                                                                       97428, 97968, 147798,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 206136, 0, 3,
                                                                       196504, 141330, 197288,
                                                                       97968, 98508, 148554,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 207144, 0, 3,
                                                                       197288, 141918, 198072,
                                                                       98508, 99048, 149310,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 208152, 0, 3,
                                                                       198072, 142506, 198856,
                                                                       99048, 99588, 150066,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 209160, 0, 3,
                                                                       198856, 143094, 199640,
                                                                       99588, 100128, 150822,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 210168, 0, 3,
                                                                       200424, 144270, 201208,
                                                                       101208, 101748, 151578,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 211176, 0, 3,
                                                                       201208, 144858, 201992,
                                                                       101748, 102288, 152334,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 212184, 0, 3,
                                                                       201992, 145446, 202776,
                                                                       102288, 102828, 153090,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 213192, 0, 3,
                                                                       202776, 146034, 203560,
                                                                       102828, 103368, 153846,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 214200, 0, 3,
                                                                       203560, 146622, 204344,
                                                                       103368, 103908, 154602,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 215208, 0, 3,
                                                                       205128, 147798, 206136,
                                                                       104988, 105663, 155358,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 216468, 0, 3,
                                                                       206136, 148554, 207144,
                                                                       105663, 106338, 156303,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 217728, 0, 3,
                                                                       207144, 149310, 208152,
                                                                       106338, 107013, 157248,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 218988, 0, 3,
                                                                       208152, 150066, 209160,
                                                                       107013, 107688, 158193,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 220248, 0, 3,
                                                                       210168, 151578, 211176,
                                                                       109038, 109713, 159138,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 221508, 0, 3,
                                                                       211176, 152334, 212184,
                                                                       109713, 110388, 160083,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 222768, 0, 3,
                                                                       212184, 153090, 213192,
                                                                       110388, 111063, 161028,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 224028, 0, 3,
                                                                       213192, 153846, 214200,
                                                                       111063, 111738, 161973,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 225288, 0, 3,
                                                                       215208, 155358, 216468,
                                                                       113088, 113913, 162918,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 226828, 0, 3,
                                                                       216468, 156303, 217728,
                                                                       113913, 114738, 164073,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 228368, 0, 3,
                                                                       217728, 157248, 218988,
                                                                       114738, 115563, 165228,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 229908, 0, 3,
                                                                       220248, 159138, 221508,
                                                                       117213, 118038, 166383,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 231448, 0, 3,
                                                                       221508, 160083, 222768,
                                                                       118038, 118863, 167538,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 232988, 0, 3,
                                                                       222768, 161028, 224028,
                                                                       118863, 119688, 168693,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 234528, 3, 121338,
                                                                       121359, 169904, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 234564, 3, 121359,
                                                                       121380, 169932, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 234600, 3, 121380,
                                                                       121401, 169960, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 234636, 3, 121401,
                                                                       121422, 169988, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 234672, 3, 121422,
                                                                       121443, 170016, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 234708, 3, 121443,
                                                                       121464, 170044, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 234744, 3, 121464,
                                                                       121485, 170072, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 234780, 3, 121485,
                                                                       121506, 170100, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 234816, 3, 121506,
                                                                       121527, 170128, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 234852, 3, 121527,
                                                                       121548, 170156, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 234888, 3, 121590,
                                                                       121611, 170240, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 234924, 3, 121611,
                                                                       121632, 170268, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 234960, 3, 121632,
                                                                       121653, 170296, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 234996, 3, 121653,
                                                                       121674, 170324, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 235032, 3, 121674,
                                                                       121695, 170352, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 235068, 3, 121695,
                                                                       121716, 170380, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 235104, 3, 121716,
                                                                       121737, 170408, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 235140, 3, 121737,
                                                                       121758, 170436, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 235176, 3, 121758,
                                                                       121779, 170464, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 235212, 3, 121779,
                                                                       121800, 170492, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 235248, 0, 3,
                                                                       234528, 169904, 234564,
                                                                       121842, 121905, 170688,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 235356, 0, 3,
                                                                       234564, 169932, 234600,
                                                                       121905, 121968, 170772,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 235464, 0, 3,
                                                                       234600, 169960, 234636,
                                                                       121968, 122031, 170856,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 235572, 0, 3,
                                                                       234636, 169988, 234672,
                                                                       122031, 122094, 170940,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 235680, 0, 3,
                                                                       234672, 170016, 234708,
                                                                       122094, 122157, 171024,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 235788, 0, 3,
                                                                       234708, 170044, 234744,
                                                                       122157, 122220, 171108,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 235896, 0, 3,
                                                                       234744, 170072, 234780,
                                                                       122220, 122283, 171192,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 236004, 0, 3,
                                                                       234780, 170100, 234816,
                                                                       122283, 122346, 171276,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 236112, 0, 3,
                                                                       234816, 170128, 234852,
                                                                       122346, 122409, 171360,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 236220, 0, 3,
                                                                       234888, 170240, 234924,
                                                                       122535, 122598, 171612,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 236328, 0, 3,
                                                                       234924, 170268, 234960,
                                                                       122598, 122661, 171696,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 236436, 0, 3,
                                                                       234960, 170296, 234996,
                                                                       122661, 122724, 171780,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 236544, 0, 3,
                                                                       234996, 170324, 235032,
                                                                       122724, 122787, 171864,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 236652, 0, 3,
                                                                       235032, 170352, 235068,
                                                                       122787, 122850, 171948,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 236760, 0, 3,
                                                                       235068, 170380, 235104,
                                                                       122850, 122913, 172032,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 236868, 0, 3,
                                                                       235104, 170408, 235140,
                                                                       122913, 122976, 172116,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 236976, 0, 3,
                                                                       235140, 170436, 235176,
                                                                       122976, 123039, 172200,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 237084, 0, 3,
                                                                       235176, 170464, 235212,
                                                                       123039, 123102, 172284,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 237192, 0, 3,
                                                                       235248, 170688, 235356,
                                                                       123228, 123354, 172704,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 237408, 0, 3,
                                                                       235356, 170772, 235464,
                                                                       123354, 123480, 172872,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 237624, 0, 3,
                                                                       235464, 170856, 235572,
                                                                       123480, 123606, 173040,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 237840, 0, 3,
                                                                       235572, 170940, 235680,
                                                                       123606, 123732, 173208,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 238056, 0, 3,
                                                                       235680, 171024, 235788,
                                                                       123732, 123858, 173376,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 238272, 0, 3,
                                                                       235788, 171108, 235896,
                                                                       123858, 123984, 173544,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 238488, 0, 3,
                                                                       235896, 171192, 236004,
                                                                       123984, 124110, 173712,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 238704, 0, 3,
                                                                       236004, 171276, 236112,
                                                                       124110, 124236, 173880,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 238920, 0, 3,
                                                                       236220, 171612, 236328,
                                                                       124488, 124614, 174384,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 239136, 0, 3,
                                                                       236328, 171696, 236436,
                                                                       124614, 124740, 174552,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 239352, 0, 3,
                                                                       236436, 171780, 236544,
                                                                       124740, 124866, 174720,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 239568, 0, 3,
                                                                       236544, 171864, 236652,
                                                                       124866, 124992, 174888,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 239784, 0, 3,
                                                                       236652, 171948, 236760,
                                                                       124992, 125118, 175056,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 240000, 0, 3,
                                                                       236760, 172032, 236868,
                                                                       125118, 125244, 175224,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 240216, 0, 3,
                                                                       236868, 172116, 236976,
                                                                       125244, 125370, 175392,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 240432, 0, 3,
                                                                       236976, 172200, 237084,
                                                                       125370, 125496, 175560,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 240648, 0, 3,
                                                                       237192, 172704, 237408,
                                                                       125748, 125958, 176288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 241008, 0, 3,
                                                                       237408, 172872, 237624,
                                                                       125958, 126168, 176568,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 241368, 0, 3,
                                                                       237624, 173040, 237840,
                                                                       126168, 126378, 176848,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 241728, 0, 3,
                                                                       237840, 173208, 238056,
                                                                       126378, 126588, 177128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 242088, 0, 3,
                                                                       238056, 173376, 238272,
                                                                       126588, 126798, 177408,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 242448, 0, 3,
                                                                       238272, 173544, 238488,
                                                                       126798, 127008, 177688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 242808, 0, 3,
                                                                       238488, 173712, 238704,
                                                                       127008, 127218, 177968,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 243168, 0, 3,
                                                                       238920, 174384, 239136,
                                                                       127638, 127848, 178808,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 243528, 0, 3,
                                                                       239136, 174552, 239352,
                                                                       127848, 128058, 179088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 243888, 0, 3,
                                                                       239352, 174720, 239568,
                                                                       128058, 128268, 179368,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 244248, 0, 3,
                                                                       239568, 174888, 239784,
                                                                       128268, 128478, 179648,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 244608, 0, 3,
                                                                       239784, 175056, 240000,
                                                                       128478, 128688, 179928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 244968, 0, 3,
                                                                       240000, 175224, 240216,
                                                                       128688, 128898, 180208,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 245328, 0, 3,
                                                                       240216, 175392, 240432,
                                                                       128898, 129108, 180488,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 245688, 0, 3,
                                                                       240648, 176288, 241008,
                                                                       129528, 129843, 181608,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 246228, 0, 3,
                                                                       241008, 176568, 241368,
                                                                       129843, 130158, 182028,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 246768, 0, 3,
                                                                       241368, 176848, 241728,
                                                                       130158, 130473, 182448,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 247308, 0, 3,
                                                                       241728, 177128, 242088,
                                                                       130473, 130788, 182868,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 247848, 0, 3,
                                                                       242088, 177408, 242448,
                                                                       130788, 131103, 183288,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 248388, 0, 3,
                                                                       242448, 177688, 242808,
                                                                       131103, 131418, 183708,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 248928, 0, 3,
                                                                       243168, 178808, 243528,
                                                                       132048, 132363, 184968,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 249468, 0, 3,
                                                                       243528, 179088, 243888,
                                                                       132363, 132678, 185388,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 250008, 0, 3,
                                                                       243888, 179368, 244248,
                                                                       132678, 132993, 185808,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 250548, 0, 3,
                                                                       244248, 179648, 244608,
                                                                       132993, 133308, 186228,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 251088, 0, 3,
                                                                       244608, 179928, 244968,
                                                                       133308, 133623, 186648,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 251628, 0, 3,
                                                                       244968, 180208, 245328,
                                                                       133623, 133938, 187068,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 252168, 0, 3,
                                                                       245688, 181608, 246228,
                                                                       134568, 135009, 188664,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 252924, 0, 3,
                                                                       246228, 182028, 246768,
                                                                       135009, 135450, 189252,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 253680, 0, 3,
                                                                       246768, 182448, 247308,
                                                                       135450, 135891, 189840,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 254436, 0, 3,
                                                                       247308, 182868, 247848,
                                                                       135891, 136332, 190428,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 255192, 0, 3,
                                                                       247848, 183288, 248388,
                                                                       136332, 136773, 191016,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 255948, 0, 3,
                                                                       248928, 184968, 249468,
                                                                       137655, 138096, 192780,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 256704, 0, 3,
                                                                       249468, 185388, 250008,
                                                                       138096, 138537, 193368,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 257460, 0, 3,
                                                                       250008, 185808, 250548,
                                                                       138537, 138978, 193956,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 258216, 0, 3,
                                                                       250548, 186228, 251088,
                                                                       138978, 139419, 194544,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 258972, 0, 3,
                                                                       251088, 186648, 251628,
                                                                       139419, 139860, 195132,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 259728, 0, 3,
                                                                       252168, 188664, 252924,
                                                                       140742, 141330, 197288,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 260736, 0, 3,
                                                                       252924, 189252, 253680,
                                                                       141330, 141918, 198072,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 261744, 0, 3,
                                                                       253680, 189840, 254436,
                                                                       141918, 142506, 198856,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 262752, 0, 3,
                                                                       254436, 190428, 255192,
                                                                       142506, 143094, 199640,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 263760, 0, 3,
                                                                       255948, 192780, 256704,
                                                                       144270, 144858, 201992,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 264768, 0, 3,
                                                                       256704, 193368, 257460,
                                                                       144858, 145446, 202776,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 265776, 0, 3,
                                                                       257460, 193956, 258216,
                                                                       145446, 146034, 203560,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 266784, 0, 3,
                                                                       258216, 194544, 258972,
                                                                       146034, 146622, 204344,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 267792, 0, 3,
                                                                       259728, 197288, 260736,
                                                                       147798, 148554, 207144,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 269088, 0, 3,
                                                                       260736, 198072, 261744,
                                                                       148554, 149310, 208152,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 270384, 0, 3,
                                                                       261744, 198856, 262752,
                                                                       149310, 150066, 209160,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 271680, 0, 3,
                                                                       263760, 201992, 264768,
                                                                       151578, 152334, 212184,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 272976, 0, 3,
                                                                       264768, 202776, 265776,
                                                                       152334, 153090, 213192,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 274272, 0, 3,
                                                                       265776, 203560, 266784,
                                                                       153090, 153846, 214200,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 275568, 0, 3,
                                                                       267792, 207144, 269088,
                                                                       155358, 156303, 217728,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 277188, 0, 3,
                                                                       269088, 208152, 270384,
                                                                       156303, 157248, 218988,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 278808, 0, 3,
                                                                       271680, 212184, 272976,
                                                                       159138, 160083, 222768,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 280428, 0, 3,
                                                                       272976, 213192, 274272,
                                                                       160083, 161028, 224028,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 282048, 0, 3,
                                                                       275568, 217728, 277188,
                                                                       162918, 164073, 228368,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 284028, 0, 3,
                                                                       278808, 222768, 280428,
                                                                       166383, 167538, 232988,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 286008, 3, 169848,
                                                                       169876, 234528, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 286053, 3, 169876,
                                                                       169904, 234564, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 286098, 3, 169904,
                                                                       169932, 234600, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 286143, 3, 169932,
                                                                       169960, 234636, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 286188, 3, 169960,
                                                                       169988, 234672, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 286233, 3, 169988,
                                                                       170016, 234708, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 286278, 3, 170016,
                                                                       170044, 234744, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 286323, 3, 170044,
                                                                       170072, 234780, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 286368, 3, 170072,
                                                                       170100, 234816, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 286413, 3, 170100,
                                                                       170128, 234852, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 286458, 3, 170184,
                                                                       170212, 234888, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 286503, 3, 170212,
                                                                       170240, 234924, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 286548, 3, 170240,
                                                                       170268, 234960, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 286593, 3, 170268,
                                                                       170296, 234996, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 286638, 3, 170296,
                                                                       170324, 235032, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 286683, 3, 170324,
                                                                       170352, 235068, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 286728, 3, 170352,
                                                                       170380, 235104, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 286773, 3, 170380,
                                                                       170408, 235140, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 286818, 3, 170408,
                                                                       170436, 235176, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 286863, 3, 170436,
                                                                       170464, 235212, ncols,
                                                                       gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 286908, 0, 3,
                                                                       286008, 234528, 286053,
                                                                       170520, 170604, 235248,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 287043, 0, 3,
                                                                       286053, 234564, 286098,
                                                                       170604, 170688, 235356,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 287178, 0, 3,
                                                                       286098, 234600, 286143,
                                                                       170688, 170772, 235464,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 287313, 0, 3,
                                                                       286143, 234636, 286188,
                                                                       170772, 170856, 235572,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 287448, 0, 3,
                                                                       286188, 234672, 286233,
                                                                       170856, 170940, 235680,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 287583, 0, 3,
                                                                       286233, 234708, 286278,
                                                                       170940, 171024, 235788,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 287718, 0, 3,
                                                                       286278, 234744, 286323,
                                                                       171024, 171108, 235896,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 287853, 0, 3,
                                                                       286323, 234780, 286368,
                                                                       171108, 171192, 236004,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 287988, 0, 3,
                                                                       286368, 234816, 286413,
                                                                       171192, 171276, 236112,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 288123, 0, 3,
                                                                       286458, 234888, 286503,
                                                                       171444, 171528, 236220,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 288258, 0, 3,
                                                                       286503, 234924, 286548,
                                                                       171528, 171612, 236328,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 288393, 0, 3,
                                                                       286548, 234960, 286593,
                                                                       171612, 171696, 236436,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 288528, 0, 3,
                                                                       286593, 234996, 286638,
                                                                       171696, 171780, 236544,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 288663, 0, 3,
                                                                       286638, 235032, 286683,
                                                                       171780, 171864, 236652,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 288798, 0, 3,
                                                                       286683, 235068, 286728,
                                                                       171864, 171948, 236760,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 288933, 0, 3,
                                                                       286728, 235104, 286773,
                                                                       171948, 172032, 236868,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 289068, 0, 3,
                                                                       286773, 235140, 286818,
                                                                       172032, 172116, 236976,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 289203, 0, 3,
                                                                       286818, 235176, 286863,
                                                                       172116, 172200, 237084,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 289338, 0, 3,
                                                                       286908, 235248, 287043,
                                                                       172368, 172536, 237192,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 289608, 0, 3,
                                                                       287043, 235356, 287178,
                                                                       172536, 172704, 237408,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 289878, 0, 3,
                                                                       287178, 235464, 287313,
                                                                       172704, 172872, 237624,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 290148, 0, 3,
                                                                       287313, 235572, 287448,
                                                                       172872, 173040, 237840,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 290418, 0, 3,
                                                                       287448, 235680, 287583,
                                                                       173040, 173208, 238056,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 290688, 0, 3,
                                                                       287583, 235788, 287718,
                                                                       173208, 173376, 238272,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 290958, 0, 3,
                                                                       287718, 235896, 287853,
                                                                       173376, 173544, 238488,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 291228, 0, 3,
                                                                       287853, 236004, 287988,
                                                                       173544, 173712, 238704,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 291498, 0, 3,
                                                                       288123, 236220, 288258,
                                                                       174048, 174216, 238920,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 291768, 0, 3,
                                                                       288258, 236328, 288393,
                                                                       174216, 174384, 239136,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 292038, 0, 3,
                                                                       288393, 236436, 288528,
                                                                       174384, 174552, 239352,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 292308, 0, 3,
                                                                       288528, 236544, 288663,
                                                                       174552, 174720, 239568,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 292578, 0, 3,
                                                                       288663, 236652, 288798,
                                                                       174720, 174888, 239784,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 292848, 0, 3,
                                                                       288798, 236760, 288933,
                                                                       174888, 175056, 240000,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 293118, 0, 3,
                                                                       288933, 236868, 289068,
                                                                       175056, 175224, 240216,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 293388, 0, 3,
                                                                       289068, 236976, 289203,
                                                                       175224, 175392, 240432,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 293658, 0, 3,
                                                                       289338, 237192, 289608,
                                                                       175728, 176008, 240648,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 294108, 0, 3,
                                                                       289608, 237408, 289878,
                                                                       176008, 176288, 241008,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 294558, 0, 3,
                                                                       289878, 237624, 290148,
                                                                       176288, 176568, 241368,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 295008, 0, 3,
                                                                       290148, 237840, 290418,
                                                                       176568, 176848, 241728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 295458, 0, 3,
                                                                       290418, 238056, 290688,
                                                                       176848, 177128, 242088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 295908, 0, 3,
                                                                       290688, 238272, 290958,
                                                                       177128, 177408, 242448,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 296358, 0, 3,
                                                                       290958, 238488, 291228,
                                                                       177408, 177688, 242808,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 296808, 0, 3,
                                                                       291498, 238920, 291768,
                                                                       178248, 178528, 243168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 297258, 0, 3,
                                                                       291768, 239136, 292038,
                                                                       178528, 178808, 243528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 297708, 0, 3,
                                                                       292038, 239352, 292308,
                                                                       178808, 179088, 243888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 298158, 0, 3,
                                                                       292308, 239568, 292578,
                                                                       179088, 179368, 244248,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 298608, 0, 3,
                                                                       292578, 239784, 292848,
                                                                       179368, 179648, 244608,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 299058, 0, 3,
                                                                       292848, 240000, 293118,
                                                                       179648, 179928, 244968,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 299508, 0, 3,
                                                                       293118, 240216, 293388,
                                                                       179928, 180208, 245328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 299958, 0, 3,
                                                                       293658, 240648, 294108,
                                                                       180768, 181188, 245688,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 300633, 0, 3,
                                                                       294108, 241008, 294558,
                                                                       181188, 181608, 246228,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 301308, 0, 3,
                                                                       294558, 241368, 295008,
                                                                       181608, 182028, 246768,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 301983, 0, 3,
                                                                       295008, 241728, 295458,
                                                                       182028, 182448, 247308,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 302658, 0, 3,
                                                                       295458, 242088, 295908,
                                                                       182448, 182868, 247848,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 303333, 0, 3,
                                                                       295908, 242448, 296358,
                                                                       182868, 183288, 248388,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 304008, 0, 3,
                                                                       296808, 243168, 297258,
                                                                       184128, 184548, 248928,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 304683, 0, 3,
                                                                       297258, 243528, 297708,
                                                                       184548, 184968, 249468,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 305358, 0, 3,
                                                                       297708, 243888, 298158,
                                                                       184968, 185388, 250008,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 306033, 0, 3,
                                                                       298158, 244248, 298608,
                                                                       185388, 185808, 250548,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 306708, 0, 3,
                                                                       298608, 244608, 299058,
                                                                       185808, 186228, 251088,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 307383, 0, 3,
                                                                       299058, 244968, 299508,
                                                                       186228, 186648, 251628,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 308058, 0, 3,
                                                                       299958, 245688, 300633,
                                                                       187488, 188076, 252168,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 309003, 0, 3,
                                                                       300633, 246228, 301308,
                                                                       188076, 188664, 252924,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 309948, 0, 3,
                                                                       301308, 246768, 301983,
                                                                       188664, 189252, 253680,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 310893, 0, 3,
                                                                       301983, 247308, 302658,
                                                                       189252, 189840, 254436,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 311838, 0, 3,
                                                                       302658, 247848, 303333,
                                                                       189840, 190428, 255192,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 312783, 0, 3,
                                                                       304008, 248928, 304683,
                                                                       191604, 192192, 255948,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 313728, 0, 3,
                                                                       304683, 249468, 305358,
                                                                       192192, 192780, 256704,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 314673, 0, 3,
                                                                       305358, 250008, 306033,
                                                                       192780, 193368, 257460,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 315618, 0, 3,
                                                                       306033, 250548, 306708,
                                                                       193368, 193956, 258216,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 316563, 0, 3,
                                                                       306708, 251088, 307383,
                                                                       193956, 194544, 258972,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 317508, 0, 3,
                                                                       308058, 252168, 309003,
                                                                       195720, 196504, 259728,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 318768, 0, 3,
                                                                       309003, 252924, 309948,
                                                                       196504, 197288, 260736,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 320028, 0, 3,
                                                                       309948, 253680, 310893,
                                                                       197288, 198072, 261744,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 321288, 0, 3,
                                                                       310893, 254436, 311838,
                                                                       198072, 198856, 262752,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 322548, 0, 3,
                                                                       312783, 255948, 313728,
                                                                       200424, 201208, 263760,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 323808, 0, 3,
                                                                       313728, 256704, 314673,
                                                                       201208, 201992, 264768,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 325068, 0, 3,
                                                                       314673, 257460, 315618,
                                                                       201992, 202776, 265776,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 326328, 0, 3,
                                                                       315618, 258216, 316563,
                                                                       202776, 203560, 266784,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 327588, 0, 3,
                                                                       317508, 259728, 318768,
                                                                       205128, 206136, 267792,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 329208, 0, 3,
                                                                       318768, 260736, 320028,
                                                                       206136, 207144, 269088,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 330828, 0, 3,
                                                                       320028, 261744, 321288,
                                                                       207144, 208152, 270384,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 332448, 0, 3,
                                                                       322548, 263760, 323808,
                                                                       210168, 211176, 271680,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 334068, 0, 3,
                                                                       323808, 264768, 325068,
                                                                       211176, 212184, 272976,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 335688, 0, 3,
                                                                       325068, 265776, 326328,
                                                                       212184, 213192, 274272,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 337308, 0, 3,
                                                                       327588, 267792, 329208,
                                                                       215208, 216468, 275568,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 339333, 0, 3,
                                                                       329208, 269088, 330828,
                                                                       216468, 217728, 277188,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 341358, 0, 3,
                                                                       332448, 271680, 334068,
                                                                       220248, 221508, 278808,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 343383, 0, 3,
                                                                       334068, 272976, 335688,
                                                                       221508, 222768, 280428,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 345408, 0, 3,
                                                                       337308, 275568, 339333,
                                                                       225288, 226828, 282048,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 347883, 0, 3,
                                                                       341358, 278808, 343383,
                                                                       229908, 231448, 284028,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 350358, 317508, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 352094, 322548, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 353830, 327588, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 356062, 332448, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 358294, 337308, 2025, ncols);

                    simdfunc::contract_primitives(buffer, 361084, 341358, 2025, ncols);

                    simdfunc::contract_primitives(buffer, 363874, 345408, 2475, ncols);

                    simdfunc::contract_primitives(buffer, 367284, 347883, 2475, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 351618, 350358, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 353354, 352094, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 355450, 353830, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 357682, 356062, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 360319, 358294, 45, 1, nmax);

        simdtrf::transform_l_inner(buffer, 363109, 361084, 45, 1, nmax);

        simdtrf::transform_l_inner(buffer, 366349, 363874, 55, 1, nmax);

        simdtrf::transform_l_inner(buffer, 369759, 367284, 55, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 370694, 351618, 355450, 17,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 372122, 353354, 357682, 17,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 373550, 355450, 360319, 17,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 375386, 357682, 363109, 17,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 377222, 360319, 366349, 17,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 379517, 363109, 369759, 17,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 381812, 370694, 373550, 17,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 384668, 372122, 375386, 17,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 387524, 373550, 377222, 17,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 391196, 375386, 379517, 17,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 394868, 381812, 387524, 17,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 399628, 384668, 391196, 17,
                                             nmax);

        simdtrf::transform_f_inner(buffer, 404388, 399628, 28, 17, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 404388, 119, nmax);

        simdtrf::transform_f_inner(buffer, 404388, 394868, 28, 17, nmax);

        simdtrf::transform_i_outer(values + 1547 * nvalues + n * npairs, nvalues, buffer, 404388,
                                   119, nmax);
    }

    for (size_t m = 0; m < 3094; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
