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


#include "SimdThreeCenterElectronRepulsionRsRecIHI.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIF.hpp"
#include "SimdTransferIG.hpp"
#include "SimdTransferIH.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKF.hpp"
#include "SimdTransferKG.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLD.hpp"
#include "SimdTransferLF.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransferMD.hpp"
#include "SimdTransferMP.hpp"
#include "SimdTransferNP.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_ihi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_ihi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 437676, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 3718 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 437676, 296304, 24242, dimensions);

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

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4958, 0, 3, 3068,
                                                                       3113, 3968, 4023, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5024, 0, 3, 3113,
                                                                       3158, 4023, 4078, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5090, 0, 3, 3158,
                                                                       3203, 4078, 4133, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5156, 0, 3, 3203,
                                                                       3248, 4133, 4188, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5222, 0, 3, 3248,
                                                                       3293, 4188, 4243, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5288, 0, 3, 3293,
                                                                       3338, 4243, 4298, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5354, 0, 3, 3338,
                                                                       3383, 4298, 4353, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5420, 0, 3, 3383,
                                                                       3428, 4353, 4408, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5486, 0, 3, 3518,
                                                                       3563, 4463, 4518, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5552, 0, 3, 3563,
                                                                       3608, 4518, 4573, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5618, 0, 3, 3608,
                                                                       3653, 4573, 4628, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5684, 0, 3, 3653,
                                                                       3698, 4628, 4683, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5750, 0, 3, 3698,
                                                                       3743, 4683, 4738, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5816, 0, 3, 3743,
                                                                       3788, 4738, 4793, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5882, 0, 3, 3788,
                                                                       3833, 4793, 4848, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5948, 0, 3, 3833,
                                                                       3878, 4848, 4903, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6014, 0, 3, 3968,
                                                                       4023, 4958, 5024, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6092, 0, 3, 4023,
                                                                       4078, 5024, 5090, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6170, 0, 3, 4078,
                                                                       4133, 5090, 5156, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6248, 0, 3, 4133,
                                                                       4188, 5156, 5222, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6326, 0, 3, 4188,
                                                                       4243, 5222, 5288, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6404, 0, 3, 4243,
                                                                       4298, 5288, 5354, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6482, 0, 3, 4298,
                                                                       4353, 5354, 5420, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6560, 0, 3, 4463,
                                                                       4518, 5486, 5552, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6638, 0, 3, 4518,
                                                                       4573, 5552, 5618, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6716, 0, 3, 4573,
                                                                       4628, 5618, 5684, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6794, 0, 3, 4628,
                                                                       4683, 5684, 5750, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6872, 0, 3, 4683,
                                                                       4738, 5750, 5816, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6950, 0, 3, 4738,
                                                                       4793, 5816, 5882, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 7028, 0, 3, 4793,
                                                                       4848, 5882, 5948, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7106, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7109, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7112, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7115, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7118, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7121, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7124, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7127, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7130, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7133, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7136, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7139, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7142, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7145, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7148, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7151, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7154, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7157, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7160, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7163, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7166, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7169, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7172, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7175, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7178, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7181, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7184, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7187, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7190, 3, 40,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7193, 3, 41,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7196, 3, 42,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7199, 3, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7202, 3, 9, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7211, 3, 10, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7220, 3, 11, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7229, 3, 12, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7238, 3, 13, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7247, 3, 14, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7256, 3, 15, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7265, 3, 16, 71,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7274, 3, 17, 74,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7283, 3, 18, 77,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7292, 3, 19, 80,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7301, 3, 20, 83,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7310, 3, 21, 86,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7319, 3, 22, 89,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7328, 3, 23, 92,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7337, 3, 28, 101,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7346, 3, 29, 104,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7355, 3, 30, 107,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7364, 3, 31, 110,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7373, 3, 32, 113,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7382, 3, 33, 116,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7391, 3, 34, 119,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7400, 3, 35, 122,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7409, 3, 36, 125,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7418, 3, 37, 128,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7427, 3, 38, 131,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7436, 3, 39, 134,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7445, 3, 40, 137,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7454, 3, 41, 140,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7463, 3, 42, 143,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7472, 3, 50, 158,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7490, 3, 53, 164,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7508, 3, 56, 170,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7526, 3, 59, 176,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7544, 3, 62, 182,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7562, 3, 65, 188,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7580, 3, 68, 194,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7598, 3, 71, 200,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7616, 3, 74, 206,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7634, 3, 77, 212,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7652, 3, 80, 218,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7670, 3, 83, 224,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7688, 3, 86, 230,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7706, 3, 89, 236,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7724, 3, 101, 254,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7742, 3, 104, 260,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7760, 3, 107, 266,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7778, 3, 110, 272,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7796, 3, 113, 278,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7814, 3, 116, 284,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7832, 3, 119, 290,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7850, 3, 122, 296,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7868, 3, 125, 302,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7886, 3, 128, 308,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7904, 3, 131, 314,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7922, 3, 134, 320,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7940, 3, 137, 326,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7958, 3, 140, 332,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 7976, 3, 158, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8006, 3, 164, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8036, 3, 170, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8066, 3, 176, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8096, 3, 182, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8126, 3, 188, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8156, 3, 194, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8186, 3, 200, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8216, 3, 206, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8246, 3, 212, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8276, 3, 218, 458,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8306, 3, 224, 468,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8336, 3, 230, 478,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8366, 3, 254, 508,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8396, 3, 260, 518,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8426, 3, 266, 528,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8456, 3, 272, 538,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8486, 3, 278, 548,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8516, 3, 284, 558,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8546, 3, 290, 568,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8576, 3, 296, 578,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8606, 3, 302, 588,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8636, 3, 308, 598,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8666, 3, 314, 608,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8696, 3, 320, 618,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8726, 3, 326, 628,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8756, 3, 358, 668,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8801, 3, 368, 683,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8846, 3, 378, 698,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8891, 3, 388, 713,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8936, 3, 398, 728,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8981, 3, 408, 743,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9026, 3, 418, 758,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9071, 3, 428, 773,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9116, 3, 438, 788,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9161, 3, 448, 803,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9206, 3, 458, 818,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9251, 3, 468, 833,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9296, 3, 508, 878,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9341, 3, 518, 893,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9386, 3, 528, 908,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9431, 3, 538, 923,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9476, 3, 548, 938,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9521, 3, 558, 953,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9566, 3, 568, 968,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9611, 3, 578, 983,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9656, 3, 588, 998,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9701, 3, 598,
                                                                       1013, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9746, 3, 608,
                                                                       1028, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9791, 3, 618,
                                                                       1043, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9836, 3, 668,
                                                                       1100, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9899, 3, 683,
                                                                       1121, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9962, 3, 698,
                                                                       1142, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10025, 3, 713,
                                                                       1163, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10088, 3, 728,
                                                                       1184, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10151, 3, 743,
                                                                       1205, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10214, 3, 758,
                                                                       1226, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10277, 3, 773,
                                                                       1247, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10340, 3, 788,
                                                                       1268, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10403, 3, 803,
                                                                       1289, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10466, 3, 818,
                                                                       1310, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10529, 3, 878,
                                                                       1373, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10592, 3, 893,
                                                                       1394, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10655, 3, 908,
                                                                       1415, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10718, 3, 923,
                                                                       1436, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10781, 3, 938,
                                                                       1457, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10844, 3, 953,
                                                                       1478, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10907, 3, 968,
                                                                       1499, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10970, 3, 983,
                                                                       1520, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11033, 3, 998,
                                                                       1541, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11096, 3, 1013,
                                                                       1562, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11159, 3, 1028,
                                                                       1583, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11222, 3, 1100,
                                                                       1660, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11306, 3, 1121,
                                                                       1688, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11390, 3, 1142,
                                                                       1716, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11474, 3, 1163,
                                                                       1744, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11558, 3, 1184,
                                                                       1772, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11642, 3, 1205,
                                                                       1800, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11726, 3, 1226,
                                                                       1828, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11810, 3, 1247,
                                                                       1856, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11894, 3, 1268,
                                                                       1884, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11978, 3, 1289,
                                                                       1912, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12062, 3, 1373,
                                                                       1996, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12146, 3, 1394,
                                                                       2024, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12230, 3, 1415,
                                                                       2052, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12314, 3, 1436,
                                                                       2080, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12398, 3, 1457,
                                                                       2108, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12482, 3, 1478,
                                                                       2136, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12566, 3, 1499,
                                                                       2164, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12650, 3, 1520,
                                                                       2192, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12734, 3, 1541,
                                                                       2220, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12818, 3, 1562,
                                                                       2248, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12902, 3, 1660,
                                                                       2348, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13010, 3, 1688,
                                                                       2384, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13118, 3, 1716,
                                                                       2420, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13226, 3, 1744,
                                                                       2456, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13334, 3, 1772,
                                                                       2492, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13442, 3, 1800,
                                                                       2528, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13550, 3, 1828,
                                                                       2564, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13658, 3, 1856,
                                                                       2600, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13766, 3, 1884,
                                                                       2636, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13874, 3, 1996,
                                                                       2744, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13982, 3, 2024,
                                                                       2780, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14090, 3, 2052,
                                                                       2816, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14198, 3, 2080,
                                                                       2852, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14306, 3, 2108,
                                                                       2888, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14414, 3, 2136,
                                                                       2924, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14522, 3, 2164,
                                                                       2960, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14630, 3, 2192,
                                                                       2996, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14738, 3, 2220,
                                                                       3032, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14846, 3, 2348,
                                                                       3158, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14981, 3, 2384,
                                                                       3203, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15116, 3, 2420,
                                                                       3248, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15251, 3, 2456,
                                                                       3293, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15386, 3, 2492,
                                                                       3338, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15521, 3, 2528,
                                                                       3383, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15656, 3, 2564,
                                                                       3428, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15791, 3, 2600,
                                                                       3473, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15926, 3, 2744,
                                                                       3608, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16061, 3, 2780,
                                                                       3653, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16196, 3, 2816,
                                                                       3698, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16331, 3, 2852,
                                                                       3743, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16466, 3, 2888,
                                                                       3788, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16601, 3, 2924,
                                                                       3833, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16736, 3, 2960,
                                                                       3878, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16871, 3, 2996,
                                                                       3923, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17006, 3, 3158,
                                                                       4078, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17171, 3, 3203,
                                                                       4133, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17336, 3, 3248,
                                                                       4188, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17501, 3, 3293,
                                                                       4243, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17666, 3, 3338,
                                                                       4298, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17831, 3, 3383,
                                                                       4353, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17996, 3, 3428,
                                                                       4408, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 18161, 3, 3608,
                                                                       4573, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 18326, 3, 3653,
                                                                       4628, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 18491, 3, 3698,
                                                                       4683, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 18656, 3, 3743,
                                                                       4738, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 18821, 3, 3788,
                                                                       4793, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 18986, 3, 3833,
                                                                       4848, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 19151, 3, 3878,
                                                                       4903, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 19316, 3, 4078,
                                                                       5090, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 19514, 3, 4133,
                                                                       5156, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 19712, 3, 4188,
                                                                       5222, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 19910, 3, 4243,
                                                                       5288, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 20108, 3, 4298,
                                                                       5354, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 20306, 3, 4353,
                                                                       5420, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 20504, 3, 4573,
                                                                       5618, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 20702, 3, 4628,
                                                                       5684, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 20900, 3, 4683,
                                                                       5750, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 21098, 3, 4738,
                                                                       5816, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 21296, 3, 4793,
                                                                       5882, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 21494, 3, 4848,
                                                                       5948, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 21692, 3, 5090,
                                                                       6170, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 21926, 3, 5156,
                                                                       6248, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 22160, 3, 5222,
                                                                       6326, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 22394, 3, 5288,
                                                                       6404, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 22628, 3, 5354,
                                                                       6482, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 22862, 3, 5618,
                                                                       6716, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 23096, 3, 5684,
                                                                       6794, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 23330, 3, 5750,
                                                                       6872, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 23564, 3, 5816,
                                                                       6950, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 23798, 3, 5882,
                                                                       7028, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24032, 3, 7, 8,
                                                                       7106, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24038, 3, 8, 9,
                                                                       7109, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24044, 3, 9, 10,
                                                                       7112, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24050, 3, 10, 11,
                                                                       7115, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24056, 3, 11, 12,
                                                                       7118, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24062, 3, 12, 13,
                                                                       7121, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24068, 3, 13, 14,
                                                                       7124, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24074, 3, 14, 15,
                                                                       7127, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24080, 3, 15, 16,
                                                                       7130, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24086, 3, 16, 17,
                                                                       7133, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24092, 3, 17, 18,
                                                                       7136, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24098, 3, 18, 19,
                                                                       7139, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24104, 3, 19, 20,
                                                                       7142, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24110, 3, 20, 21,
                                                                       7145, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24116, 3, 21, 22,
                                                                       7148, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24122, 3, 22, 23,
                                                                       7151, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24128, 3, 26, 27,
                                                                       7154, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24134, 3, 27, 28,
                                                                       7157, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24140, 3, 28, 29,
                                                                       7160, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24146, 3, 29, 30,
                                                                       7163, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24152, 3, 30, 31,
                                                                       7166, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24158, 3, 31, 32,
                                                                       7169, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24164, 3, 32, 33,
                                                                       7172, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24170, 3, 33, 34,
                                                                       7175, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24176, 3, 34, 35,
                                                                       7178, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24182, 3, 35, 36,
                                                                       7181, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24188, 3, 36, 37,
                                                                       7184, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24194, 3, 37, 38,
                                                                       7187, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24200, 3, 38, 39,
                                                                       7190, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24206, 3, 39, 40,
                                                                       7193, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24212, 3, 40, 41,
                                                                       7196, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 24218, 3, 41, 42,
                                                                       7199, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24224, 0, 3,
                                                                       24032, 7106, 24038, 7202,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24242, 0, 3,
                                                                       24038, 7109, 24044, 7211,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24260, 0, 3,
                                                                       24044, 7112, 24050, 7220,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24278, 0, 3,
                                                                       24050, 7115, 24056, 7229,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24296, 0, 3,
                                                                       24056, 7118, 24062, 7238,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24314, 0, 3,
                                                                       24062, 7121, 24068, 7247,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24332, 0, 3,
                                                                       24068, 7124, 24074, 7256,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24350, 0, 3,
                                                                       24074, 7127, 24080, 7265,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24368, 0, 3,
                                                                       24080, 7130, 24086, 7274,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24386, 0, 3,
                                                                       24086, 7133, 24092, 7283,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24404, 0, 3,
                                                                       24092, 7136, 24098, 7292,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24422, 0, 3,
                                                                       24098, 7139, 24104, 7301,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24440, 0, 3,
                                                                       24104, 7142, 24110, 7310,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24458, 0, 3,
                                                                       24110, 7145, 24116, 7319,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24476, 0, 3,
                                                                       24116, 7148, 24122, 7328,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24494, 0, 3,
                                                                       24128, 7154, 24134, 7337,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24512, 0, 3,
                                                                       24134, 7157, 24140, 7346,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24530, 0, 3,
                                                                       24140, 7160, 24146, 7355,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24548, 0, 3,
                                                                       24146, 7163, 24152, 7364,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24566, 0, 3,
                                                                       24152, 7166, 24158, 7373,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24584, 0, 3,
                                                                       24158, 7169, 24164, 7382,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24602, 0, 3,
                                                                       24164, 7172, 24170, 7391,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24620, 0, 3,
                                                                       24170, 7175, 24176, 7400,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24638, 0, 3,
                                                                       24176, 7178, 24182, 7409,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24656, 0, 3,
                                                                       24182, 7181, 24188, 7418,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24674, 0, 3,
                                                                       24188, 7184, 24194, 7427,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24692, 0, 3,
                                                                       24194, 7187, 24200, 7436,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24710, 0, 3,
                                                                       24200, 7190, 24206, 7445,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24728, 0, 3,
                                                                       24206, 7193, 24212, 7454,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 24746, 0, 3,
                                                                       24212, 7196, 24218, 7463,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24764, 0, 3,
                                                                       24224, 7202, 24242, 146,
                                                                       152, 7472, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24800, 0, 3,
                                                                       24242, 7211, 24260, 152,
                                                                       158, 7490, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24836, 0, 3,
                                                                       24260, 7220, 24278, 158,
                                                                       164, 7508, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24872, 0, 3,
                                                                       24278, 7229, 24296, 164,
                                                                       170, 7526, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24908, 0, 3,
                                                                       24296, 7238, 24314, 170,
                                                                       176, 7544, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24944, 0, 3,
                                                                       24314, 7247, 24332, 176,
                                                                       182, 7562, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24980, 0, 3,
                                                                       24332, 7256, 24350, 182,
                                                                       188, 7580, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25016, 0, 3,
                                                                       24350, 7265, 24368, 188,
                                                                       194, 7598, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25052, 0, 3,
                                                                       24368, 7274, 24386, 194,
                                                                       200, 7616, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25088, 0, 3,
                                                                       24386, 7283, 24404, 200,
                                                                       206, 7634, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25124, 0, 3,
                                                                       24404, 7292, 24422, 206,
                                                                       212, 7652, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25160, 0, 3,
                                                                       24422, 7301, 24440, 212,
                                                                       218, 7670, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25196, 0, 3,
                                                                       24440, 7310, 24458, 218,
                                                                       224, 7688, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25232, 0, 3,
                                                                       24458, 7319, 24476, 224,
                                                                       230, 7706, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25268, 0, 3,
                                                                       24494, 7337, 24512, 242,
                                                                       248, 7724, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25304, 0, 3,
                                                                       24512, 7346, 24530, 248,
                                                                       254, 7742, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25340, 0, 3,
                                                                       24530, 7355, 24548, 254,
                                                                       260, 7760, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25376, 0, 3,
                                                                       24548, 7364, 24566, 260,
                                                                       266, 7778, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25412, 0, 3,
                                                                       24566, 7373, 24584, 266,
                                                                       272, 7796, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25448, 0, 3,
                                                                       24584, 7382, 24602, 272,
                                                                       278, 7814, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25484, 0, 3,
                                                                       24602, 7391, 24620, 278,
                                                                       284, 7832, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25520, 0, 3,
                                                                       24620, 7400, 24638, 284,
                                                                       290, 7850, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25556, 0, 3,
                                                                       24638, 7409, 24656, 290,
                                                                       296, 7868, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25592, 0, 3,
                                                                       24656, 7418, 24674, 296,
                                                                       302, 7886, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25628, 0, 3,
                                                                       24674, 7427, 24692, 302,
                                                                       308, 7904, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25664, 0, 3,
                                                                       24692, 7436, 24710, 308,
                                                                       314, 7922, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25700, 0, 3,
                                                                       24710, 7445, 24728, 314,
                                                                       320, 7940, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 25736, 0, 3,
                                                                       24728, 7454, 24746, 320,
                                                                       326, 7958, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25772, 0, 3,
                                                                       24764, 7472, 24800, 338,
                                                                       348, 7976, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25832, 0, 3,
                                                                       24800, 7490, 24836, 348,
                                                                       358, 8006, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25892, 0, 3,
                                                                       24836, 7508, 24872, 358,
                                                                       368, 8036, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25952, 0, 3,
                                                                       24872, 7526, 24908, 368,
                                                                       378, 8066, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26012, 0, 3,
                                                                       24908, 7544, 24944, 378,
                                                                       388, 8096, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26072, 0, 3,
                                                                       24944, 7562, 24980, 388,
                                                                       398, 8126, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26132, 0, 3,
                                                                       24980, 7580, 25016, 398,
                                                                       408, 8156, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26192, 0, 3,
                                                                       25016, 7598, 25052, 408,
                                                                       418, 8186, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26252, 0, 3,
                                                                       25052, 7616, 25088, 418,
                                                                       428, 8216, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26312, 0, 3,
                                                                       25088, 7634, 25124, 428,
                                                                       438, 8246, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26372, 0, 3,
                                                                       25124, 7652, 25160, 438,
                                                                       448, 8276, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26432, 0, 3,
                                                                       25160, 7670, 25196, 448,
                                                                       458, 8306, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26492, 0, 3,
                                                                       25196, 7688, 25232, 458,
                                                                       468, 8336, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26552, 0, 3,
                                                                       25268, 7724, 25304, 488,
                                                                       498, 8366, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26612, 0, 3,
                                                                       25304, 7742, 25340, 498,
                                                                       508, 8396, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26672, 0, 3,
                                                                       25340, 7760, 25376, 508,
                                                                       518, 8426, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26732, 0, 3,
                                                                       25376, 7778, 25412, 518,
                                                                       528, 8456, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26792, 0, 3,
                                                                       25412, 7796, 25448, 528,
                                                                       538, 8486, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26852, 0, 3,
                                                                       25448, 7814, 25484, 538,
                                                                       548, 8516, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26912, 0, 3,
                                                                       25484, 7832, 25520, 548,
                                                                       558, 8546, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 26972, 0, 3,
                                                                       25520, 7850, 25556, 558,
                                                                       568, 8576, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 27032, 0, 3,
                                                                       25556, 7868, 25592, 568,
                                                                       578, 8606, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 27092, 0, 3,
                                                                       25592, 7886, 25628, 578,
                                                                       588, 8636, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 27152, 0, 3,
                                                                       25628, 7904, 25664, 588,
                                                                       598, 8666, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 27212, 0, 3,
                                                                       25664, 7922, 25700, 598,
                                                                       608, 8696, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 27272, 0, 3,
                                                                       25700, 7940, 25736, 608,
                                                                       618, 8726, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27332, 0, 3,
                                                                       25772, 7976, 25832, 638,
                                                                       653, 8756, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27422, 0, 3,
                                                                       25832, 8006, 25892, 653,
                                                                       668, 8801, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27512, 0, 3,
                                                                       25892, 8036, 25952, 668,
                                                                       683, 8846, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27602, 0, 3,
                                                                       25952, 8066, 26012, 683,
                                                                       698, 8891, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27692, 0, 3,
                                                                       26012, 8096, 26072, 698,
                                                                       713, 8936, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27782, 0, 3,
                                                                       26072, 8126, 26132, 713,
                                                                       728, 8981, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27872, 0, 3,
                                                                       26132, 8156, 26192, 728,
                                                                       743, 9026, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27962, 0, 3,
                                                                       26192, 8186, 26252, 743,
                                                                       758, 9071, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 28052, 0, 3,
                                                                       26252, 8216, 26312, 758,
                                                                       773, 9116, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 28142, 0, 3,
                                                                       26312, 8246, 26372, 773,
                                                                       788, 9161, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 28232, 0, 3,
                                                                       26372, 8276, 26432, 788,
                                                                       803, 9206, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 28322, 0, 3,
                                                                       26432, 8306, 26492, 803,
                                                                       818, 9251, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 28412, 0, 3,
                                                                       26552, 8366, 26612, 848,
                                                                       863, 9296, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 28502, 0, 3,
                                                                       26612, 8396, 26672, 863,
                                                                       878, 9341, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 28592, 0, 3,
                                                                       26672, 8426, 26732, 878,
                                                                       893, 9386, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 28682, 0, 3,
                                                                       26732, 8456, 26792, 893,
                                                                       908, 9431, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 28772, 0, 3,
                                                                       26792, 8486, 26852, 908,
                                                                       923, 9476, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 28862, 0, 3,
                                                                       26852, 8516, 26912, 923,
                                                                       938, 9521, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 28952, 0, 3,
                                                                       26912, 8546, 26972, 938,
                                                                       953, 9566, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 29042, 0, 3,
                                                                       26972, 8576, 27032, 953,
                                                                       968, 9611, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 29132, 0, 3,
                                                                       27032, 8606, 27092, 968,
                                                                       983, 9656, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 29222, 0, 3,
                                                                       27092, 8636, 27152, 983,
                                                                       998, 9701, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 29312, 0, 3,
                                                                       27152, 8666, 27212, 998,
                                                                       1013, 9746, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 29402, 0, 3,
                                                                       27212, 8696, 27272, 1013,
                                                                       1028, 9791, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 29492, 0, 3,
                                                                       27332, 8756, 27422, 1058,
                                                                       1079, 9836, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 29618, 0, 3,
                                                                       27422, 8801, 27512, 1079,
                                                                       1100, 9899, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 29744, 0, 3,
                                                                       27512, 8846, 27602, 1100,
                                                                       1121, 9962, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 29870, 0, 3,
                                                                       27602, 8891, 27692, 1121,
                                                                       1142, 10025, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 29996, 0, 3,
                                                                       27692, 8936, 27782, 1142,
                                                                       1163, 10088, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 30122, 0, 3,
                                                                       27782, 8981, 27872, 1163,
                                                                       1184, 10151, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 30248, 0, 3,
                                                                       27872, 9026, 27962, 1184,
                                                                       1205, 10214, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 30374, 0, 3,
                                                                       27962, 9071, 28052, 1205,
                                                                       1226, 10277, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 30500, 0, 3,
                                                                       28052, 9116, 28142, 1226,
                                                                       1247, 10340, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 30626, 0, 3,
                                                                       28142, 9161, 28232, 1247,
                                                                       1268, 10403, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 30752, 0, 3,
                                                                       28232, 9206, 28322, 1268,
                                                                       1289, 10466, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 30878, 0, 3,
                                                                       28412, 9296, 28502, 1331,
                                                                       1352, 10529, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 31004, 0, 3,
                                                                       28502, 9341, 28592, 1352,
                                                                       1373, 10592, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 31130, 0, 3,
                                                                       28592, 9386, 28682, 1373,
                                                                       1394, 10655, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 31256, 0, 3,
                                                                       28682, 9431, 28772, 1394,
                                                                       1415, 10718, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 31382, 0, 3,
                                                                       28772, 9476, 28862, 1415,
                                                                       1436, 10781, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 31508, 0, 3,
                                                                       28862, 9521, 28952, 1436,
                                                                       1457, 10844, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 31634, 0, 3,
                                                                       28952, 9566, 29042, 1457,
                                                                       1478, 10907, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 31760, 0, 3,
                                                                       29042, 9611, 29132, 1478,
                                                                       1499, 10970, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 31886, 0, 3,
                                                                       29132, 9656, 29222, 1499,
                                                                       1520, 11033, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 32012, 0, 3,
                                                                       29222, 9701, 29312, 1520,
                                                                       1541, 11096, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 32138, 0, 3,
                                                                       29312, 9746, 29402, 1541,
                                                                       1562, 11159, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 32264, 0, 3,
                                                                       29492, 9836, 29618, 1604,
                                                                       1632, 11222, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 32432, 0, 3,
                                                                       29618, 9899, 29744, 1632,
                                                                       1660, 11306, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 32600, 0, 3,
                                                                       29744, 9962, 29870, 1660,
                                                                       1688, 11390, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 32768, 0, 3,
                                                                       29870, 10025, 29996, 1688,
                                                                       1716, 11474, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 32936, 0, 3,
                                                                       29996, 10088, 30122, 1716,
                                                                       1744, 11558, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 33104, 0, 3,
                                                                       30122, 10151, 30248, 1744,
                                                                       1772, 11642, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 33272, 0, 3,
                                                                       30248, 10214, 30374, 1772,
                                                                       1800, 11726, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 33440, 0, 3,
                                                                       30374, 10277, 30500, 1800,
                                                                       1828, 11810, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 33608, 0, 3,
                                                                       30500, 10340, 30626, 1828,
                                                                       1856, 11894, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 33776, 0, 3,
                                                                       30626, 10403, 30752, 1856,
                                                                       1884, 11978, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 33944, 0, 3,
                                                                       30878, 10529, 31004, 1940,
                                                                       1968, 12062, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 34112, 0, 3,
                                                                       31004, 10592, 31130, 1968,
                                                                       1996, 12146, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 34280, 0, 3,
                                                                       31130, 10655, 31256, 1996,
                                                                       2024, 12230, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 34448, 0, 3,
                                                                       31256, 10718, 31382, 2024,
                                                                       2052, 12314, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 34616, 0, 3,
                                                                       31382, 10781, 31508, 2052,
                                                                       2080, 12398, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 34784, 0, 3,
                                                                       31508, 10844, 31634, 2080,
                                                                       2108, 12482, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 34952, 0, 3,
                                                                       31634, 10907, 31760, 2108,
                                                                       2136, 12566, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 35120, 0, 3,
                                                                       31760, 10970, 31886, 2136,
                                                                       2164, 12650, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 35288, 0, 3,
                                                                       31886, 11033, 32012, 2164,
                                                                       2192, 12734, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 35456, 0, 3,
                                                                       32012, 11096, 32138, 2192,
                                                                       2220, 12818, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 35624, 0, 3,
                                                                       32264, 11222, 32432, 2276,
                                                                       2312, 12902, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 35840, 0, 3,
                                                                       32432, 11306, 32600, 2312,
                                                                       2348, 13010, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 36056, 0, 3,
                                                                       32600, 11390, 32768, 2348,
                                                                       2384, 13118, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 36272, 0, 3,
                                                                       32768, 11474, 32936, 2384,
                                                                       2420, 13226, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 36488, 0, 3,
                                                                       32936, 11558, 33104, 2420,
                                                                       2456, 13334, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 36704, 0, 3,
                                                                       33104, 11642, 33272, 2456,
                                                                       2492, 13442, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 36920, 0, 3,
                                                                       33272, 11726, 33440, 2492,
                                                                       2528, 13550, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 37136, 0, 3,
                                                                       33440, 11810, 33608, 2528,
                                                                       2564, 13658, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 37352, 0, 3,
                                                                       33608, 11894, 33776, 2564,
                                                                       2600, 13766, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 37568, 0, 3,
                                                                       33944, 12062, 34112, 2672,
                                                                       2708, 13874, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 37784, 0, 3,
                                                                       34112, 12146, 34280, 2708,
                                                                       2744, 13982, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 38000, 0, 3,
                                                                       34280, 12230, 34448, 2744,
                                                                       2780, 14090, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 38216, 0, 3,
                                                                       34448, 12314, 34616, 2780,
                                                                       2816, 14198, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 38432, 0, 3,
                                                                       34616, 12398, 34784, 2816,
                                                                       2852, 14306, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 38648, 0, 3,
                                                                       34784, 12482, 34952, 2852,
                                                                       2888, 14414, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 38864, 0, 3,
                                                                       34952, 12566, 35120, 2888,
                                                                       2924, 14522, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 39080, 0, 3,
                                                                       35120, 12650, 35288, 2924,
                                                                       2960, 14630, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 39296, 0, 3,
                                                                       35288, 12734, 35456, 2960,
                                                                       2996, 14738, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 39512, 0, 3,
                                                                       35624, 12902, 35840, 3068,
                                                                       3113, 14846, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 39782, 0, 3,
                                                                       35840, 13010, 36056, 3113,
                                                                       3158, 14981, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 40052, 0, 3,
                                                                       36056, 13118, 36272, 3158,
                                                                       3203, 15116, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 40322, 0, 3,
                                                                       36272, 13226, 36488, 3203,
                                                                       3248, 15251, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 40592, 0, 3,
                                                                       36488, 13334, 36704, 3248,
                                                                       3293, 15386, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 40862, 0, 3,
                                                                       36704, 13442, 36920, 3293,
                                                                       3338, 15521, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 41132, 0, 3,
                                                                       36920, 13550, 37136, 3338,
                                                                       3383, 15656, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 41402, 0, 3,
                                                                       37136, 13658, 37352, 3383,
                                                                       3428, 15791, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 41672, 0, 3,
                                                                       37568, 13874, 37784, 3518,
                                                                       3563, 15926, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 41942, 0, 3,
                                                                       37784, 13982, 38000, 3563,
                                                                       3608, 16061, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 42212, 0, 3,
                                                                       38000, 14090, 38216, 3608,
                                                                       3653, 16196, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 42482, 0, 3,
                                                                       38216, 14198, 38432, 3653,
                                                                       3698, 16331, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 42752, 0, 3,
                                                                       38432, 14306, 38648, 3698,
                                                                       3743, 16466, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 43022, 0, 3,
                                                                       38648, 14414, 38864, 3743,
                                                                       3788, 16601, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 43292, 0, 3,
                                                                       38864, 14522, 39080, 3788,
                                                                       3833, 16736, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 43562, 0, 3,
                                                                       39080, 14630, 39296, 3833,
                                                                       3878, 16871, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 43832, 0, 3,
                                                                       39512, 14846, 39782, 3968,
                                                                       4023, 17006, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 44162, 0, 3,
                                                                       39782, 14981, 40052, 4023,
                                                                       4078, 17171, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 44492, 0, 3,
                                                                       40052, 15116, 40322, 4078,
                                                                       4133, 17336, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 44822, 0, 3,
                                                                       40322, 15251, 40592, 4133,
                                                                       4188, 17501, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 45152, 0, 3,
                                                                       40592, 15386, 40862, 4188,
                                                                       4243, 17666, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 45482, 0, 3,
                                                                       40862, 15521, 41132, 4243,
                                                                       4298, 17831, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 45812, 0, 3,
                                                                       41132, 15656, 41402, 4298,
                                                                       4353, 17996, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 46142, 0, 3,
                                                                       41672, 15926, 41942, 4463,
                                                                       4518, 18161, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 46472, 0, 3,
                                                                       41942, 16061, 42212, 4518,
                                                                       4573, 18326, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 46802, 0, 3,
                                                                       42212, 16196, 42482, 4573,
                                                                       4628, 18491, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 47132, 0, 3,
                                                                       42482, 16331, 42752, 4628,
                                                                       4683, 18656, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 47462, 0, 3,
                                                                       42752, 16466, 43022, 4683,
                                                                       4738, 18821, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 47792, 0, 3,
                                                                       43022, 16601, 43292, 4738,
                                                                       4793, 18986, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 48122, 0, 3,
                                                                       43292, 16736, 43562, 4793,
                                                                       4848, 19151, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 48452, 0, 3,
                                                                       43832, 17006, 44162, 4958,
                                                                       5024, 19316, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 48848, 0, 3,
                                                                       44162, 17171, 44492, 5024,
                                                                       5090, 19514, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 49244, 0, 3,
                                                                       44492, 17336, 44822, 5090,
                                                                       5156, 19712, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 49640, 0, 3,
                                                                       44822, 17501, 45152, 5156,
                                                                       5222, 19910, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 50036, 0, 3,
                                                                       45152, 17666, 45482, 5222,
                                                                       5288, 20108, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 50432, 0, 3,
                                                                       45482, 17831, 45812, 5288,
                                                                       5354, 20306, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 50828, 0, 3,
                                                                       46142, 18161, 46472, 5486,
                                                                       5552, 20504, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 51224, 0, 3,
                                                                       46472, 18326, 46802, 5552,
                                                                       5618, 20702, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 51620, 0, 3,
                                                                       46802, 18491, 47132, 5618,
                                                                       5684, 20900, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 52016, 0, 3,
                                                                       47132, 18656, 47462, 5684,
                                                                       5750, 21098, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 52412, 0, 3,
                                                                       47462, 18821, 47792, 5750,
                                                                       5816, 21296, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 52808, 0, 3,
                                                                       47792, 18986, 48122, 5816,
                                                                       5882, 21494, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 53204, 0, 3,
                                                                       48452, 19316, 48848, 6014,
                                                                       6092, 21692, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 53672, 0, 3,
                                                                       48848, 19514, 49244, 6092,
                                                                       6170, 21926, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 54140, 0, 3,
                                                                       49244, 19712, 49640, 6170,
                                                                       6248, 22160, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 54608, 0, 3,
                                                                       49640, 19910, 50036, 6248,
                                                                       6326, 22394, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 55076, 0, 3,
                                                                       50036, 20108, 50432, 6326,
                                                                       6404, 22628, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 55544, 0, 3,
                                                                       50828, 20504, 51224, 6560,
                                                                       6638, 22862, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 56012, 0, 3,
                                                                       51224, 20702, 51620, 6638,
                                                                       6716, 23096, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 56480, 0, 3,
                                                                       51620, 20900, 52016, 6716,
                                                                       6794, 23330, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 56948, 0, 3,
                                                                       52016, 21098, 52412, 6794,
                                                                       6872, 23564, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 57416, 0, 3,
                                                                       52412, 21296, 52808, 6872,
                                                                       6950, 23798, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 57884, 3, 7106,
                                                                       7109, 24044, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 57894, 3, 7109,
                                                                       7112, 24050, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 57904, 3, 7112,
                                                                       7115, 24056, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 57914, 3, 7115,
                                                                       7118, 24062, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 57924, 3, 7118,
                                                                       7121, 24068, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 57934, 3, 7121,
                                                                       7124, 24074, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 57944, 3, 7124,
                                                                       7127, 24080, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 57954, 3, 7127,
                                                                       7130, 24086, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 57964, 3, 7130,
                                                                       7133, 24092, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 57974, 3, 7133,
                                                                       7136, 24098, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 57984, 3, 7136,
                                                                       7139, 24104, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 57994, 3, 7139,
                                                                       7142, 24110, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 58004, 3, 7142,
                                                                       7145, 24116, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 58014, 3, 7145,
                                                                       7148, 24122, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 58024, 3, 7154,
                                                                       7157, 24140, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 58034, 3, 7157,
                                                                       7160, 24146, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 58044, 3, 7160,
                                                                       7163, 24152, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 58054, 3, 7163,
                                                                       7166, 24158, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 58064, 3, 7166,
                                                                       7169, 24164, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 58074, 3, 7169,
                                                                       7172, 24170, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 58084, 3, 7172,
                                                                       7175, 24176, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 58094, 3, 7175,
                                                                       7178, 24182, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 58104, 3, 7178,
                                                                       7181, 24188, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 58114, 3, 7181,
                                                                       7184, 24194, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 58124, 3, 7184,
                                                                       7187, 24200, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 58134, 3, 7187,
                                                                       7190, 24206, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 58144, 3, 7190,
                                                                       7193, 24212, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 58154, 3, 7193,
                                                                       7196, 24218, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58164, 0, 3,
                                                                       57884, 24044, 57894, 7202,
                                                                       7211, 24260, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58194, 0, 3,
                                                                       57894, 24050, 57904, 7211,
                                                                       7220, 24278, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58224, 0, 3,
                                                                       57904, 24056, 57914, 7220,
                                                                       7229, 24296, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58254, 0, 3,
                                                                       57914, 24062, 57924, 7229,
                                                                       7238, 24314, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58284, 0, 3,
                                                                       57924, 24068, 57934, 7238,
                                                                       7247, 24332, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58314, 0, 3,
                                                                       57934, 24074, 57944, 7247,
                                                                       7256, 24350, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58344, 0, 3,
                                                                       57944, 24080, 57954, 7256,
                                                                       7265, 24368, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58374, 0, 3,
                                                                       57954, 24086, 57964, 7265,
                                                                       7274, 24386, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58404, 0, 3,
                                                                       57964, 24092, 57974, 7274,
                                                                       7283, 24404, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58434, 0, 3,
                                                                       57974, 24098, 57984, 7283,
                                                                       7292, 24422, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58464, 0, 3,
                                                                       57984, 24104, 57994, 7292,
                                                                       7301, 24440, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58494, 0, 3,
                                                                       57994, 24110, 58004, 7301,
                                                                       7310, 24458, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58524, 0, 3,
                                                                       58004, 24116, 58014, 7310,
                                                                       7319, 24476, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58554, 0, 3,
                                                                       58024, 24140, 58034, 7337,
                                                                       7346, 24530, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58584, 0, 3,
                                                                       58034, 24146, 58044, 7346,
                                                                       7355, 24548, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58614, 0, 3,
                                                                       58044, 24152, 58054, 7355,
                                                                       7364, 24566, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58644, 0, 3,
                                                                       58054, 24158, 58064, 7364,
                                                                       7373, 24584, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58674, 0, 3,
                                                                       58064, 24164, 58074, 7373,
                                                                       7382, 24602, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58704, 0, 3,
                                                                       58074, 24170, 58084, 7382,
                                                                       7391, 24620, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58734, 0, 3,
                                                                       58084, 24176, 58094, 7391,
                                                                       7400, 24638, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58764, 0, 3,
                                                                       58094, 24182, 58104, 7400,
                                                                       7409, 24656, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58794, 0, 3,
                                                                       58104, 24188, 58114, 7409,
                                                                       7418, 24674, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58824, 0, 3,
                                                                       58114, 24194, 58124, 7418,
                                                                       7427, 24692, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58854, 0, 3,
                                                                       58124, 24200, 58134, 7427,
                                                                       7436, 24710, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58884, 0, 3,
                                                                       58134, 24206, 58144, 7436,
                                                                       7445, 24728, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 58914, 0, 3,
                                                                       58144, 24212, 58154, 7445,
                                                                       7454, 24746, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 58944, 0, 3,
                                                                       58164, 24260, 58194, 7472,
                                                                       7490, 24836, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 59004, 0, 3,
                                                                       58194, 24278, 58224, 7490,
                                                                       7508, 24872, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 59064, 0, 3,
                                                                       58224, 24296, 58254, 7508,
                                                                       7526, 24908, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 59124, 0, 3,
                                                                       58254, 24314, 58284, 7526,
                                                                       7544, 24944, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 59184, 0, 3,
                                                                       58284, 24332, 58314, 7544,
                                                                       7562, 24980, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 59244, 0, 3,
                                                                       58314, 24350, 58344, 7562,
                                                                       7580, 25016, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 59304, 0, 3,
                                                                       58344, 24368, 58374, 7580,
                                                                       7598, 25052, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 59364, 0, 3,
                                                                       58374, 24386, 58404, 7598,
                                                                       7616, 25088, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 59424, 0, 3,
                                                                       58404, 24404, 58434, 7616,
                                                                       7634, 25124, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 59484, 0, 3,
                                                                       58434, 24422, 58464, 7634,
                                                                       7652, 25160, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 59544, 0, 3,
                                                                       58464, 24440, 58494, 7652,
                                                                       7670, 25196, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 59604, 0, 3,
                                                                       58494, 24458, 58524, 7670,
                                                                       7688, 25232, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 59664, 0, 3,
                                                                       58554, 24530, 58584, 7724,
                                                                       7742, 25340, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 59724, 0, 3,
                                                                       58584, 24548, 58614, 7742,
                                                                       7760, 25376, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 59784, 0, 3,
                                                                       58614, 24566, 58644, 7760,
                                                                       7778, 25412, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 59844, 0, 3,
                                                                       58644, 24584, 58674, 7778,
                                                                       7796, 25448, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 59904, 0, 3,
                                                                       58674, 24602, 58704, 7796,
                                                                       7814, 25484, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 59964, 0, 3,
                                                                       58704, 24620, 58734, 7814,
                                                                       7832, 25520, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 60024, 0, 3,
                                                                       58734, 24638, 58764, 7832,
                                                                       7850, 25556, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 60084, 0, 3,
                                                                       58764, 24656, 58794, 7850,
                                                                       7868, 25592, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 60144, 0, 3,
                                                                       58794, 24674, 58824, 7868,
                                                                       7886, 25628, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 60204, 0, 3,
                                                                       58824, 24692, 58854, 7886,
                                                                       7904, 25664, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 60264, 0, 3,
                                                                       58854, 24710, 58884, 7904,
                                                                       7922, 25700, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 60324, 0, 3,
                                                                       58884, 24728, 58914, 7922,
                                                                       7940, 25736, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 60384, 0, 3,
                                                                       58944, 24836, 59004, 7976,
                                                                       8006, 25892, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 60484, 0, 3,
                                                                       59004, 24872, 59064, 8006,
                                                                       8036, 25952, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 60584, 0, 3,
                                                                       59064, 24908, 59124, 8036,
                                                                       8066, 26012, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 60684, 0, 3,
                                                                       59124, 24944, 59184, 8066,
                                                                       8096, 26072, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 60784, 0, 3,
                                                                       59184, 24980, 59244, 8096,
                                                                       8126, 26132, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 60884, 0, 3,
                                                                       59244, 25016, 59304, 8126,
                                                                       8156, 26192, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 60984, 0, 3,
                                                                       59304, 25052, 59364, 8156,
                                                                       8186, 26252, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 61084, 0, 3,
                                                                       59364, 25088, 59424, 8186,
                                                                       8216, 26312, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 61184, 0, 3,
                                                                       59424, 25124, 59484, 8216,
                                                                       8246, 26372, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 61284, 0, 3,
                                                                       59484, 25160, 59544, 8246,
                                                                       8276, 26432, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 61384, 0, 3,
                                                                       59544, 25196, 59604, 8276,
                                                                       8306, 26492, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 61484, 0, 3,
                                                                       59664, 25340, 59724, 8366,
                                                                       8396, 26672, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 61584, 0, 3,
                                                                       59724, 25376, 59784, 8396,
                                                                       8426, 26732, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 61684, 0, 3,
                                                                       59784, 25412, 59844, 8426,
                                                                       8456, 26792, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 61784, 0, 3,
                                                                       59844, 25448, 59904, 8456,
                                                                       8486, 26852, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 61884, 0, 3,
                                                                       59904, 25484, 59964, 8486,
                                                                       8516, 26912, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 61984, 0, 3,
                                                                       59964, 25520, 60024, 8516,
                                                                       8546, 26972, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 62084, 0, 3,
                                                                       60024, 25556, 60084, 8546,
                                                                       8576, 27032, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 62184, 0, 3,
                                                                       60084, 25592, 60144, 8576,
                                                                       8606, 27092, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 62284, 0, 3,
                                                                       60144, 25628, 60204, 8606,
                                                                       8636, 27152, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 62384, 0, 3,
                                                                       60204, 25664, 60264, 8636,
                                                                       8666, 27212, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 62484, 0, 3,
                                                                       60264, 25700, 60324, 8666,
                                                                       8696, 27272, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 62584, 0, 3,
                                                                       60384, 25892, 60484, 8756,
                                                                       8801, 27512, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 62734, 0, 3,
                                                                       60484, 25952, 60584, 8801,
                                                                       8846, 27602, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 62884, 0, 3,
                                                                       60584, 26012, 60684, 8846,
                                                                       8891, 27692, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 63034, 0, 3,
                                                                       60684, 26072, 60784, 8891,
                                                                       8936, 27782, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 63184, 0, 3,
                                                                       60784, 26132, 60884, 8936,
                                                                       8981, 27872, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 63334, 0, 3,
                                                                       60884, 26192, 60984, 8981,
                                                                       9026, 27962, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 63484, 0, 3,
                                                                       60984, 26252, 61084, 9026,
                                                                       9071, 28052, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 63634, 0, 3,
                                                                       61084, 26312, 61184, 9071,
                                                                       9116, 28142, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 63784, 0, 3,
                                                                       61184, 26372, 61284, 9116,
                                                                       9161, 28232, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 63934, 0, 3,
                                                                       61284, 26432, 61384, 9161,
                                                                       9206, 28322, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 64084, 0, 3,
                                                                       61484, 26672, 61584, 9296,
                                                                       9341, 28592, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 64234, 0, 3,
                                                                       61584, 26732, 61684, 9341,
                                                                       9386, 28682, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 64384, 0, 3,
                                                                       61684, 26792, 61784, 9386,
                                                                       9431, 28772, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 64534, 0, 3,
                                                                       61784, 26852, 61884, 9431,
                                                                       9476, 28862, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 64684, 0, 3,
                                                                       61884, 26912, 61984, 9476,
                                                                       9521, 28952, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 64834, 0, 3,
                                                                       61984, 26972, 62084, 9521,
                                                                       9566, 29042, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 64984, 0, 3,
                                                                       62084, 27032, 62184, 9566,
                                                                       9611, 29132, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 65134, 0, 3,
                                                                       62184, 27092, 62284, 9611,
                                                                       9656, 29222, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 65284, 0, 3,
                                                                       62284, 27152, 62384, 9656,
                                                                       9701, 29312, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 65434, 0, 3,
                                                                       62384, 27212, 62484, 9701,
                                                                       9746, 29402, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 65584, 0, 3,
                                                                       62584, 27512, 62734, 9836,
                                                                       9899, 29744, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 65794, 0, 3,
                                                                       62734, 27602, 62884, 9899,
                                                                       9962, 29870, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 66004, 0, 3,
                                                                       62884, 27692, 63034, 9962,
                                                                       10025, 29996, ncols,
                                                                       gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 66214, 0, 3,
                                                                       63034, 27782, 63184,
                                                                       10025, 10088, 30122,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 66424, 0, 3,
                                                                       63184, 27872, 63334,
                                                                       10088, 10151, 30248,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 66634, 0, 3,
                                                                       63334, 27962, 63484,
                                                                       10151, 10214, 30374,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 66844, 0, 3,
                                                                       63484, 28052, 63634,
                                                                       10214, 10277, 30500,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 67054, 0, 3,
                                                                       63634, 28142, 63784,
                                                                       10277, 10340, 30626,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 67264, 0, 3,
                                                                       63784, 28232, 63934,
                                                                       10340, 10403, 30752,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 67474, 0, 3,
                                                                       64084, 28592, 64234,
                                                                       10529, 10592, 31130,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 67684, 0, 3,
                                                                       64234, 28682, 64384,
                                                                       10592, 10655, 31256,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 67894, 0, 3,
                                                                       64384, 28772, 64534,
                                                                       10655, 10718, 31382,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 68104, 0, 3,
                                                                       64534, 28862, 64684,
                                                                       10718, 10781, 31508,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 68314, 0, 3,
                                                                       64684, 28952, 64834,
                                                                       10781, 10844, 31634,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 68524, 0, 3,
                                                                       64834, 29042, 64984,
                                                                       10844, 10907, 31760,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 68734, 0, 3,
                                                                       64984, 29132, 65134,
                                                                       10907, 10970, 31886,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 68944, 0, 3,
                                                                       65134, 29222, 65284,
                                                                       10970, 11033, 32012,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 69154, 0, 3,
                                                                       65284, 29312, 65434,
                                                                       11033, 11096, 32138,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 69364, 0, 3,
                                                                       65584, 29744, 65794,
                                                                       11222, 11306, 32600,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 69644, 0, 3,
                                                                       65794, 29870, 66004,
                                                                       11306, 11390, 32768,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 69924, 0, 3,
                                                                       66004, 29996, 66214,
                                                                       11390, 11474, 32936,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 70204, 0, 3,
                                                                       66214, 30122, 66424,
                                                                       11474, 11558, 33104,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 70484, 0, 3,
                                                                       66424, 30248, 66634,
                                                                       11558, 11642, 33272,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 70764, 0, 3,
                                                                       66634, 30374, 66844,
                                                                       11642, 11726, 33440,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 71044, 0, 3,
                                                                       66844, 30500, 67054,
                                                                       11726, 11810, 33608,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 71324, 0, 3,
                                                                       67054, 30626, 67264,
                                                                       11810, 11894, 33776,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 71604, 0, 3,
                                                                       67474, 31130, 67684,
                                                                       12062, 12146, 34280,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 71884, 0, 3,
                                                                       67684, 31256, 67894,
                                                                       12146, 12230, 34448,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 72164, 0, 3,
                                                                       67894, 31382, 68104,
                                                                       12230, 12314, 34616,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 72444, 0, 3,
                                                                       68104, 31508, 68314,
                                                                       12314, 12398, 34784,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 72724, 0, 3,
                                                                       68314, 31634, 68524,
                                                                       12398, 12482, 34952,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 73004, 0, 3,
                                                                       68524, 31760, 68734,
                                                                       12482, 12566, 35120,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 73284, 0, 3,
                                                                       68734, 31886, 68944,
                                                                       12566, 12650, 35288,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 73564, 0, 3,
                                                                       68944, 32012, 69154,
                                                                       12650, 12734, 35456,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 73844, 0, 3,
                                                                       69364, 32600, 69644,
                                                                       12902, 13010, 36056,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 74204, 0, 3,
                                                                       69644, 32768, 69924,
                                                                       13010, 13118, 36272,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 74564, 0, 3,
                                                                       69924, 32936, 70204,
                                                                       13118, 13226, 36488,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 74924, 0, 3,
                                                                       70204, 33104, 70484,
                                                                       13226, 13334, 36704,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 75284, 0, 3,
                                                                       70484, 33272, 70764,
                                                                       13334, 13442, 36920,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 75644, 0, 3,
                                                                       70764, 33440, 71044,
                                                                       13442, 13550, 37136,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 76004, 0, 3,
                                                                       71044, 33608, 71324,
                                                                       13550, 13658, 37352,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 76364, 0, 3,
                                                                       71604, 34280, 71884,
                                                                       13874, 13982, 38000,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 76724, 0, 3,
                                                                       71884, 34448, 72164,
                                                                       13982, 14090, 38216,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 77084, 0, 3,
                                                                       72164, 34616, 72444,
                                                                       14090, 14198, 38432,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 77444, 0, 3,
                                                                       72444, 34784, 72724,
                                                                       14198, 14306, 38648,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 77804, 0, 3,
                                                                       72724, 34952, 73004,
                                                                       14306, 14414, 38864,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 78164, 0, 3,
                                                                       73004, 35120, 73284,
                                                                       14414, 14522, 39080,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 78524, 0, 3,
                                                                       73284, 35288, 73564,
                                                                       14522, 14630, 39296,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 78884, 0, 3,
                                                                       73844, 36056, 74204,
                                                                       14846, 14981, 40052,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 79334, 0, 3,
                                                                       74204, 36272, 74564,
                                                                       14981, 15116, 40322,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 79784, 0, 3,
                                                                       74564, 36488, 74924,
                                                                       15116, 15251, 40592,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 80234, 0, 3,
                                                                       74924, 36704, 75284,
                                                                       15251, 15386, 40862,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 80684, 0, 3,
                                                                       75284, 36920, 75644,
                                                                       15386, 15521, 41132,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 81134, 0, 3,
                                                                       75644, 37136, 76004,
                                                                       15521, 15656, 41402,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 81584, 0, 3,
                                                                       76364, 38000, 76724,
                                                                       15926, 16061, 42212,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 82034, 0, 3,
                                                                       76724, 38216, 77084,
                                                                       16061, 16196, 42482,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 82484, 0, 3,
                                                                       77084, 38432, 77444,
                                                                       16196, 16331, 42752,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 82934, 0, 3,
                                                                       77444, 38648, 77804,
                                                                       16331, 16466, 43022,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 83384, 0, 3,
                                                                       77804, 38864, 78164,
                                                                       16466, 16601, 43292,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 83834, 0, 3,
                                                                       78164, 39080, 78524,
                                                                       16601, 16736, 43562,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 84284, 0, 3,
                                                                       78884, 40052, 79334,
                                                                       17006, 17171, 44492,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 84834, 0, 3,
                                                                       79334, 40322, 79784,
                                                                       17171, 17336, 44822,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 85384, 0, 3,
                                                                       79784, 40592, 80234,
                                                                       17336, 17501, 45152,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 85934, 0, 3,
                                                                       80234, 40862, 80684,
                                                                       17501, 17666, 45482,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 86484, 0, 3,
                                                                       80684, 41132, 81134,
                                                                       17666, 17831, 45812,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 87034, 0, 3,
                                                                       81584, 42212, 82034,
                                                                       18161, 18326, 46802,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 87584, 0, 3,
                                                                       82034, 42482, 82484,
                                                                       18326, 18491, 47132,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 88134, 0, 3,
                                                                       82484, 42752, 82934,
                                                                       18491, 18656, 47462,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 88684, 0, 3,
                                                                       82934, 43022, 83384,
                                                                       18656, 18821, 47792,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 89234, 0, 3,
                                                                       83384, 43292, 83834,
                                                                       18821, 18986, 48122,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 89784, 0, 3,
                                                                       84284, 44492, 84834,
                                                                       19316, 19514, 49244,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 90444, 0, 3,
                                                                       84834, 44822, 85384,
                                                                       19514, 19712, 49640,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 91104, 0, 3,
                                                                       85384, 45152, 85934,
                                                                       19712, 19910, 50036,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 91764, 0, 3,
                                                                       85934, 45482, 86484,
                                                                       19910, 20108, 50432,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 92424, 0, 3,
                                                                       87034, 46802, 87584,
                                                                       20504, 20702, 51620,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 93084, 0, 3,
                                                                       87584, 47132, 88134,
                                                                       20702, 20900, 52016,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 93744, 0, 3,
                                                                       88134, 47462, 88684,
                                                                       20900, 21098, 52412,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 94404, 0, 3,
                                                                       88684, 47792, 89234,
                                                                       21098, 21296, 52808,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 95064, 0, 3,
                                                                       89784, 49244, 90444,
                                                                       21692, 21926, 54140,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 95844, 0, 3,
                                                                       90444, 49640, 91104,
                                                                       21926, 22160, 54608,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 96624, 0, 3,
                                                                       91104, 50036, 91764,
                                                                       22160, 22394, 55076,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 97404, 0, 3,
                                                                       92424, 51620, 93084,
                                                                       22862, 23096, 56480,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 98184, 0, 3,
                                                                       93084, 52016, 93744,
                                                                       23096, 23330, 56948,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 98964, 0, 3,
                                                                       93744, 52412, 94404,
                                                                       23330, 23564, 57416,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 99744, 3, 24032,
                                                                       24038, 57884, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 99759, 3, 24038,
                                                                       24044, 57894, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 99774, 3, 24044,
                                                                       24050, 57904, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 99789, 3, 24050,
                                                                       24056, 57914, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 99804, 3, 24056,
                                                                       24062, 57924, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 99819, 3, 24062,
                                                                       24068, 57934, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 99834, 3, 24068,
                                                                       24074, 57944, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 99849, 3, 24074,
                                                                       24080, 57954, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 99864, 3, 24080,
                                                                       24086, 57964, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 99879, 3, 24086,
                                                                       24092, 57974, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 99894, 3, 24092,
                                                                       24098, 57984, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 99909, 3, 24098,
                                                                       24104, 57994, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 99924, 3, 24104,
                                                                       24110, 58004, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 99939, 3, 24110,
                                                                       24116, 58014, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 99954, 3, 24128,
                                                                       24134, 58024, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 99969, 3, 24134,
                                                                       24140, 58034, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 99984, 3, 24140,
                                                                       24146, 58044, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 99999, 3, 24146,
                                                                       24152, 58054, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 100014, 3, 24152,
                                                                       24158, 58064, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 100029, 3, 24158,
                                                                       24164, 58074, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 100044, 3, 24164,
                                                                       24170, 58084, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 100059, 3, 24170,
                                                                       24176, 58094, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 100074, 3, 24176,
                                                                       24182, 58104, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 100089, 3, 24182,
                                                                       24188, 58114, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 100104, 3, 24188,
                                                                       24194, 58124, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 100119, 3, 24194,
                                                                       24200, 58134, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 100134, 3, 24200,
                                                                       24206, 58144, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 100149, 3, 24206,
                                                                       24212, 58154, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 100164, 0, 3,
                                                                       99744, 57884, 99759,
                                                                       24224, 24242, 58164,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 100209, 0, 3,
                                                                       99759, 57894, 99774,
                                                                       24242, 24260, 58194,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 100254, 0, 3,
                                                                       99774, 57904, 99789,
                                                                       24260, 24278, 58224,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 100299, 0, 3,
                                                                       99789, 57914, 99804,
                                                                       24278, 24296, 58254,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 100344, 0, 3,
                                                                       99804, 57924, 99819,
                                                                       24296, 24314, 58284,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 100389, 0, 3,
                                                                       99819, 57934, 99834,
                                                                       24314, 24332, 58314,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 100434, 0, 3,
                                                                       99834, 57944, 99849,
                                                                       24332, 24350, 58344,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 100479, 0, 3,
                                                                       99849, 57954, 99864,
                                                                       24350, 24368, 58374,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 100524, 0, 3,
                                                                       99864, 57964, 99879,
                                                                       24368, 24386, 58404,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 100569, 0, 3,
                                                                       99879, 57974, 99894,
                                                                       24386, 24404, 58434,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 100614, 0, 3,
                                                                       99894, 57984, 99909,
                                                                       24404, 24422, 58464,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 100659, 0, 3,
                                                                       99909, 57994, 99924,
                                                                       24422, 24440, 58494,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 100704, 0, 3,
                                                                       99924, 58004, 99939,
                                                                       24440, 24458, 58524,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 100749, 0, 3,
                                                                       99954, 58024, 99969,
                                                                       24494, 24512, 58554,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 100794, 0, 3,
                                                                       99969, 58034, 99984,
                                                                       24512, 24530, 58584,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 100839, 0, 3,
                                                                       99984, 58044, 99999,
                                                                       24530, 24548, 58614,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 100884, 0, 3,
                                                                       99999, 58054, 100014,
                                                                       24548, 24566, 58644,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 100929, 0, 3,
                                                                       100014, 58064, 100029,
                                                                       24566, 24584, 58674,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 100974, 0, 3,
                                                                       100029, 58074, 100044,
                                                                       24584, 24602, 58704,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 101019, 0, 3,
                                                                       100044, 58084, 100059,
                                                                       24602, 24620, 58734,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 101064, 0, 3,
                                                                       100059, 58094, 100074,
                                                                       24620, 24638, 58764,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 101109, 0, 3,
                                                                       100074, 58104, 100089,
                                                                       24638, 24656, 58794,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 101154, 0, 3,
                                                                       100089, 58114, 100104,
                                                                       24656, 24674, 58824,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 101199, 0, 3,
                                                                       100104, 58124, 100119,
                                                                       24674, 24692, 58854,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 101244, 0, 3,
                                                                       100119, 58134, 100134,
                                                                       24692, 24710, 58884,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 101289, 0, 3,
                                                                       100134, 58144, 100149,
                                                                       24710, 24728, 58914,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 101334, 0, 3,
                                                                       100164, 58164, 100209,
                                                                       24764, 24800, 58944,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 101424, 0, 3,
                                                                       100209, 58194, 100254,
                                                                       24800, 24836, 59004,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 101514, 0, 3,
                                                                       100254, 58224, 100299,
                                                                       24836, 24872, 59064,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 101604, 0, 3,
                                                                       100299, 58254, 100344,
                                                                       24872, 24908, 59124,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 101694, 0, 3,
                                                                       100344, 58284, 100389,
                                                                       24908, 24944, 59184,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 101784, 0, 3,
                                                                       100389, 58314, 100434,
                                                                       24944, 24980, 59244,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 101874, 0, 3,
                                                                       100434, 58344, 100479,
                                                                       24980, 25016, 59304,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 101964, 0, 3,
                                                                       100479, 58374, 100524,
                                                                       25016, 25052, 59364,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 102054, 0, 3,
                                                                       100524, 58404, 100569,
                                                                       25052, 25088, 59424,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 102144, 0, 3,
                                                                       100569, 58434, 100614,
                                                                       25088, 25124, 59484,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 102234, 0, 3,
                                                                       100614, 58464, 100659,
                                                                       25124, 25160, 59544,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 102324, 0, 3,
                                                                       100659, 58494, 100704,
                                                                       25160, 25196, 59604,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 102414, 0, 3,
                                                                       100749, 58554, 100794,
                                                                       25268, 25304, 59664,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 102504, 0, 3,
                                                                       100794, 58584, 100839,
                                                                       25304, 25340, 59724,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 102594, 0, 3,
                                                                       100839, 58614, 100884,
                                                                       25340, 25376, 59784,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 102684, 0, 3,
                                                                       100884, 58644, 100929,
                                                                       25376, 25412, 59844,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 102774, 0, 3,
                                                                       100929, 58674, 100974,
                                                                       25412, 25448, 59904,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 102864, 0, 3,
                                                                       100974, 58704, 101019,
                                                                       25448, 25484, 59964,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 102954, 0, 3,
                                                                       101019, 58734, 101064,
                                                                       25484, 25520, 60024,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 103044, 0, 3,
                                                                       101064, 58764, 101109,
                                                                       25520, 25556, 60084,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 103134, 0, 3,
                                                                       101109, 58794, 101154,
                                                                       25556, 25592, 60144,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 103224, 0, 3,
                                                                       101154, 58824, 101199,
                                                                       25592, 25628, 60204,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 103314, 0, 3,
                                                                       101199, 58854, 101244,
                                                                       25628, 25664, 60264,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 103404, 0, 3,
                                                                       101244, 58884, 101289,
                                                                       25664, 25700, 60324,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 103494, 0, 3,
                                                                       101334, 58944, 101424,
                                                                       25772, 25832, 60384,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 103644, 0, 3,
                                                                       101424, 59004, 101514,
                                                                       25832, 25892, 60484,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 103794, 0, 3,
                                                                       101514, 59064, 101604,
                                                                       25892, 25952, 60584,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 103944, 0, 3,
                                                                       101604, 59124, 101694,
                                                                       25952, 26012, 60684,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 104094, 0, 3,
                                                                       101694, 59184, 101784,
                                                                       26012, 26072, 60784,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 104244, 0, 3,
                                                                       101784, 59244, 101874,
                                                                       26072, 26132, 60884,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 104394, 0, 3,
                                                                       101874, 59304, 101964,
                                                                       26132, 26192, 60984,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 104544, 0, 3,
                                                                       101964, 59364, 102054,
                                                                       26192, 26252, 61084,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 104694, 0, 3,
                                                                       102054, 59424, 102144,
                                                                       26252, 26312, 61184,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 104844, 0, 3,
                                                                       102144, 59484, 102234,
                                                                       26312, 26372, 61284,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 104994, 0, 3,
                                                                       102234, 59544, 102324,
                                                                       26372, 26432, 61384,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 105144, 0, 3,
                                                                       102414, 59664, 102504,
                                                                       26552, 26612, 61484,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 105294, 0, 3,
                                                                       102504, 59724, 102594,
                                                                       26612, 26672, 61584,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 105444, 0, 3,
                                                                       102594, 59784, 102684,
                                                                       26672, 26732, 61684,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 105594, 0, 3,
                                                                       102684, 59844, 102774,
                                                                       26732, 26792, 61784,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 105744, 0, 3,
                                                                       102774, 59904, 102864,
                                                                       26792, 26852, 61884,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 105894, 0, 3,
                                                                       102864, 59964, 102954,
                                                                       26852, 26912, 61984,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 106044, 0, 3,
                                                                       102954, 60024, 103044,
                                                                       26912, 26972, 62084,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 106194, 0, 3,
                                                                       103044, 60084, 103134,
                                                                       26972, 27032, 62184,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 106344, 0, 3,
                                                                       103134, 60144, 103224,
                                                                       27032, 27092, 62284,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 106494, 0, 3,
                                                                       103224, 60204, 103314,
                                                                       27092, 27152, 62384,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 106644, 0, 3,
                                                                       103314, 60264, 103404,
                                                                       27152, 27212, 62484,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 106794, 0, 3,
                                                                       103494, 60384, 103644,
                                                                       27332, 27422, 62584,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 107019, 0, 3,
                                                                       103644, 60484, 103794,
                                                                       27422, 27512, 62734,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 107244, 0, 3,
                                                                       103794, 60584, 103944,
                                                                       27512, 27602, 62884,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 107469, 0, 3,
                                                                       103944, 60684, 104094,
                                                                       27602, 27692, 63034,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 107694, 0, 3,
                                                                       104094, 60784, 104244,
                                                                       27692, 27782, 63184,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 107919, 0, 3,
                                                                       104244, 60884, 104394,
                                                                       27782, 27872, 63334,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 108144, 0, 3,
                                                                       104394, 60984, 104544,
                                                                       27872, 27962, 63484,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 108369, 0, 3,
                                                                       104544, 61084, 104694,
                                                                       27962, 28052, 63634,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 108594, 0, 3,
                                                                       104694, 61184, 104844,
                                                                       28052, 28142, 63784,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 108819, 0, 3,
                                                                       104844, 61284, 104994,
                                                                       28142, 28232, 63934,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 109044, 0, 3,
                                                                       105144, 61484, 105294,
                                                                       28412, 28502, 64084,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 109269, 0, 3,
                                                                       105294, 61584, 105444,
                                                                       28502, 28592, 64234,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 109494, 0, 3,
                                                                       105444, 61684, 105594,
                                                                       28592, 28682, 64384,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 109719, 0, 3,
                                                                       105594, 61784, 105744,
                                                                       28682, 28772, 64534,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 109944, 0, 3,
                                                                       105744, 61884, 105894,
                                                                       28772, 28862, 64684,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 110169, 0, 3,
                                                                       105894, 61984, 106044,
                                                                       28862, 28952, 64834,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 110394, 0, 3,
                                                                       106044, 62084, 106194,
                                                                       28952, 29042, 64984,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 110619, 0, 3,
                                                                       106194, 62184, 106344,
                                                                       29042, 29132, 65134,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 110844, 0, 3,
                                                                       106344, 62284, 106494,
                                                                       29132, 29222, 65284,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 111069, 0, 3,
                                                                       106494, 62384, 106644,
                                                                       29222, 29312, 65434,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 111294, 0, 3,
                                                                       106794, 62584, 107019,
                                                                       29492, 29618, 65584,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 111609, 0, 3,
                                                                       107019, 62734, 107244,
                                                                       29618, 29744, 65794,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 111924, 0, 3,
                                                                       107244, 62884, 107469,
                                                                       29744, 29870, 66004,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 112239, 0, 3,
                                                                       107469, 63034, 107694,
                                                                       29870, 29996, 66214,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 112554, 0, 3,
                                                                       107694, 63184, 107919,
                                                                       29996, 30122, 66424,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 112869, 0, 3,
                                                                       107919, 63334, 108144,
                                                                       30122, 30248, 66634,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 113184, 0, 3,
                                                                       108144, 63484, 108369,
                                                                       30248, 30374, 66844,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 113499, 0, 3,
                                                                       108369, 63634, 108594,
                                                                       30374, 30500, 67054,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 113814, 0, 3,
                                                                       108594, 63784, 108819,
                                                                       30500, 30626, 67264,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 114129, 0, 3,
                                                                       109044, 64084, 109269,
                                                                       30878, 31004, 67474,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 114444, 0, 3,
                                                                       109269, 64234, 109494,
                                                                       31004, 31130, 67684,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 114759, 0, 3,
                                                                       109494, 64384, 109719,
                                                                       31130, 31256, 67894,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 115074, 0, 3,
                                                                       109719, 64534, 109944,
                                                                       31256, 31382, 68104,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 115389, 0, 3,
                                                                       109944, 64684, 110169,
                                                                       31382, 31508, 68314,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 115704, 0, 3,
                                                                       110169, 64834, 110394,
                                                                       31508, 31634, 68524,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 116019, 0, 3,
                                                                       110394, 64984, 110619,
                                                                       31634, 31760, 68734,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 116334, 0, 3,
                                                                       110619, 65134, 110844,
                                                                       31760, 31886, 68944,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 116649, 0, 3,
                                                                       110844, 65284, 111069,
                                                                       31886, 32012, 69154,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 116964, 0, 3,
                                                                       111294, 65584, 111609,
                                                                       32264, 32432, 69364,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 117384, 0, 3,
                                                                       111609, 65794, 111924,
                                                                       32432, 32600, 69644,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 117804, 0, 3,
                                                                       111924, 66004, 112239,
                                                                       32600, 32768, 69924,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 118224, 0, 3,
                                                                       112239, 66214, 112554,
                                                                       32768, 32936, 70204,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 118644, 0, 3,
                                                                       112554, 66424, 112869,
                                                                       32936, 33104, 70484,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 119064, 0, 3,
                                                                       112869, 66634, 113184,
                                                                       33104, 33272, 70764,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 119484, 0, 3,
                                                                       113184, 66844, 113499,
                                                                       33272, 33440, 71044,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 119904, 0, 3,
                                                                       113499, 67054, 113814,
                                                                       33440, 33608, 71324,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 120324, 0, 3,
                                                                       114129, 67474, 114444,
                                                                       33944, 34112, 71604,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 120744, 0, 3,
                                                                       114444, 67684, 114759,
                                                                       34112, 34280, 71884,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 121164, 0, 3,
                                                                       114759, 67894, 115074,
                                                                       34280, 34448, 72164,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 121584, 0, 3,
                                                                       115074, 68104, 115389,
                                                                       34448, 34616, 72444,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 122004, 0, 3,
                                                                       115389, 68314, 115704,
                                                                       34616, 34784, 72724,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 122424, 0, 3,
                                                                       115704, 68524, 116019,
                                                                       34784, 34952, 73004,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 122844, 0, 3,
                                                                       116019, 68734, 116334,
                                                                       34952, 35120, 73284,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 123264, 0, 3,
                                                                       116334, 68944, 116649,
                                                                       35120, 35288, 73564,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 123684, 0, 3,
                                                                       116964, 69364, 117384,
                                                                       35624, 35840, 73844,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 124224, 0, 3,
                                                                       117384, 69644, 117804,
                                                                       35840, 36056, 74204,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 124764, 0, 3,
                                                                       117804, 69924, 118224,
                                                                       36056, 36272, 74564,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 125304, 0, 3,
                                                                       118224, 70204, 118644,
                                                                       36272, 36488, 74924,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 125844, 0, 3,
                                                                       118644, 70484, 119064,
                                                                       36488, 36704, 75284,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 126384, 0, 3,
                                                                       119064, 70764, 119484,
                                                                       36704, 36920, 75644,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 126924, 0, 3,
                                                                       119484, 71044, 119904,
                                                                       36920, 37136, 76004,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 127464, 0, 3,
                                                                       120324, 71604, 120744,
                                                                       37568, 37784, 76364,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 128004, 0, 3,
                                                                       120744, 71884, 121164,
                                                                       37784, 38000, 76724,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 128544, 0, 3,
                                                                       121164, 72164, 121584,
                                                                       38000, 38216, 77084,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 129084, 0, 3,
                                                                       121584, 72444, 122004,
                                                                       38216, 38432, 77444,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 129624, 0, 3,
                                                                       122004, 72724, 122424,
                                                                       38432, 38648, 77804,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 130164, 0, 3,
                                                                       122424, 73004, 122844,
                                                                       38648, 38864, 78164,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 130704, 0, 3,
                                                                       122844, 73284, 123264,
                                                                       38864, 39080, 78524,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 131244, 0, 3,
                                                                       123684, 73844, 124224,
                                                                       39512, 39782, 78884,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 131919, 0, 3,
                                                                       124224, 74204, 124764,
                                                                       39782, 40052, 79334,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 132594, 0, 3,
                                                                       124764, 74564, 125304,
                                                                       40052, 40322, 79784,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 133269, 0, 3,
                                                                       125304, 74924, 125844,
                                                                       40322, 40592, 80234,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 133944, 0, 3,
                                                                       125844, 75284, 126384,
                                                                       40592, 40862, 80684,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 134619, 0, 3,
                                                                       126384, 75644, 126924,
                                                                       40862, 41132, 81134,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 135294, 0, 3,
                                                                       127464, 76364, 128004,
                                                                       41672, 41942, 81584,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 135969, 0, 3,
                                                                       128004, 76724, 128544,
                                                                       41942, 42212, 82034,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 136644, 0, 3,
                                                                       128544, 77084, 129084,
                                                                       42212, 42482, 82484,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 137319, 0, 3,
                                                                       129084, 77444, 129624,
                                                                       42482, 42752, 82934,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 137994, 0, 3,
                                                                       129624, 77804, 130164,
                                                                       42752, 43022, 83384,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 138669, 0, 3,
                                                                       130164, 78164, 130704,
                                                                       43022, 43292, 83834,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 139344, 0, 3,
                                                                       131244, 78884, 131919,
                                                                       43832, 44162, 84284,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 140169, 0, 3,
                                                                       131919, 79334, 132594,
                                                                       44162, 44492, 84834,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 140994, 0, 3,
                                                                       132594, 79784, 133269,
                                                                       44492, 44822, 85384,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 141819, 0, 3,
                                                                       133269, 80234, 133944,
                                                                       44822, 45152, 85934,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 142644, 0, 3,
                                                                       133944, 80684, 134619,
                                                                       45152, 45482, 86484,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 143469, 0, 3,
                                                                       135294, 81584, 135969,
                                                                       46142, 46472, 87034,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 144294, 0, 3,
                                                                       135969, 82034, 136644,
                                                                       46472, 46802, 87584,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 145119, 0, 3,
                                                                       136644, 82484, 137319,
                                                                       46802, 47132, 88134,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 145944, 0, 3,
                                                                       137319, 82934, 137994,
                                                                       47132, 47462, 88684,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 146769, 0, 3,
                                                                       137994, 83384, 138669,
                                                                       47462, 47792, 89234,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 147594, 0, 3,
                                                                       139344, 84284, 140169,
                                                                       48452, 48848, 89784,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 148584, 0, 3,
                                                                       140169, 84834, 140994,
                                                                       48848, 49244, 90444,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 149574, 0, 3,
                                                                       140994, 85384, 141819,
                                                                       49244, 49640, 91104,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 150564, 0, 3,
                                                                       141819, 85934, 142644,
                                                                       49640, 50036, 91764,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 151554, 0, 3,
                                                                       143469, 87034, 144294,
                                                                       50828, 51224, 92424,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 152544, 0, 3,
                                                                       144294, 87584, 145119,
                                                                       51224, 51620, 93084,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 153534, 0, 3,
                                                                       145119, 88134, 145944,
                                                                       51620, 52016, 93744,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 154524, 0, 3,
                                                                       145944, 88684, 146769,
                                                                       52016, 52412, 94404,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 155514, 0, 3,
                                                                       147594, 89784, 148584,
                                                                       53204, 53672, 95064,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 156684, 0, 3,
                                                                       148584, 90444, 149574,
                                                                       53672, 54140, 95844,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 157854, 0, 3,
                                                                       149574, 91104, 150564,
                                                                       54140, 54608, 96624,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 159024, 0, 3,
                                                                       151554, 92424, 152544,
                                                                       55544, 56012, 97404,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 160194, 0, 3,
                                                                       152544, 93084, 153534,
                                                                       56012, 56480, 98184,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 161364, 0, 3,
                                                                       153534, 93744, 154524,
                                                                       56480, 56948, 98964,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162534, 3, 57884,
                                                                       57894, 99774, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162555, 3, 57894,
                                                                       57904, 99789, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162576, 3, 57904,
                                                                       57914, 99804, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162597, 3, 57914,
                                                                       57924, 99819, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162618, 3, 57924,
                                                                       57934, 99834, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162639, 3, 57934,
                                                                       57944, 99849, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162660, 3, 57944,
                                                                       57954, 99864, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162681, 3, 57954,
                                                                       57964, 99879, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162702, 3, 57964,
                                                                       57974, 99894, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162723, 3, 57974,
                                                                       57984, 99909, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162744, 3, 57984,
                                                                       57994, 99924, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162765, 3, 57994,
                                                                       58004, 99939, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162786, 3, 58024,
                                                                       58034, 99984, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162807, 3, 58034,
                                                                       58044, 99999, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162828, 3, 58044,
                                                                       58054, 100014, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162849, 3, 58054,
                                                                       58064, 100029, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162870, 3, 58064,
                                                                       58074, 100044, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162891, 3, 58074,
                                                                       58084, 100059, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162912, 3, 58084,
                                                                       58094, 100074, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162933, 3, 58094,
                                                                       58104, 100089, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162954, 3, 58104,
                                                                       58114, 100104, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162975, 3, 58114,
                                                                       58124, 100119, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 162996, 3, 58124,
                                                                       58134, 100134, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 163017, 3, 58134,
                                                                       58144, 100149, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 163038, 0, 3,
                                                                       162534, 99774, 162555,
                                                                       58164, 58194, 100254,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 163101, 0, 3,
                                                                       162555, 99789, 162576,
                                                                       58194, 58224, 100299,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 163164, 0, 3,
                                                                       162576, 99804, 162597,
                                                                       58224, 58254, 100344,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 163227, 0, 3,
                                                                       162597, 99819, 162618,
                                                                       58254, 58284, 100389,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 163290, 0, 3,
                                                                       162618, 99834, 162639,
                                                                       58284, 58314, 100434,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 163353, 0, 3,
                                                                       162639, 99849, 162660,
                                                                       58314, 58344, 100479,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 163416, 0, 3,
                                                                       162660, 99864, 162681,
                                                                       58344, 58374, 100524,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 163479, 0, 3,
                                                                       162681, 99879, 162702,
                                                                       58374, 58404, 100569,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 163542, 0, 3,
                                                                       162702, 99894, 162723,
                                                                       58404, 58434, 100614,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 163605, 0, 3,
                                                                       162723, 99909, 162744,
                                                                       58434, 58464, 100659,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 163668, 0, 3,
                                                                       162744, 99924, 162765,
                                                                       58464, 58494, 100704,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 163731, 0, 3,
                                                                       162786, 99984, 162807,
                                                                       58554, 58584, 100839,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 163794, 0, 3,
                                                                       162807, 99999, 162828,
                                                                       58584, 58614, 100884,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 163857, 0, 3,
                                                                       162828, 100014, 162849,
                                                                       58614, 58644, 100929,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 163920, 0, 3,
                                                                       162849, 100029, 162870,
                                                                       58644, 58674, 100974,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 163983, 0, 3,
                                                                       162870, 100044, 162891,
                                                                       58674, 58704, 101019,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 164046, 0, 3,
                                                                       162891, 100059, 162912,
                                                                       58704, 58734, 101064,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 164109, 0, 3,
                                                                       162912, 100074, 162933,
                                                                       58734, 58764, 101109,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 164172, 0, 3,
                                                                       162933, 100089, 162954,
                                                                       58764, 58794, 101154,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 164235, 0, 3,
                                                                       162954, 100104, 162975,
                                                                       58794, 58824, 101199,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 164298, 0, 3,
                                                                       162975, 100119, 162996,
                                                                       58824, 58854, 101244,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 164361, 0, 3,
                                                                       162996, 100134, 163017,
                                                                       58854, 58884, 101289,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 164424, 0, 3,
                                                                       163038, 100254, 163101,
                                                                       58944, 59004, 101514,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 164550, 0, 3,
                                                                       163101, 100299, 163164,
                                                                       59004, 59064, 101604,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 164676, 0, 3,
                                                                       163164, 100344, 163227,
                                                                       59064, 59124, 101694,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 164802, 0, 3,
                                                                       163227, 100389, 163290,
                                                                       59124, 59184, 101784,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 164928, 0, 3,
                                                                       163290, 100434, 163353,
                                                                       59184, 59244, 101874,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 165054, 0, 3,
                                                                       163353, 100479, 163416,
                                                                       59244, 59304, 101964,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 165180, 0, 3,
                                                                       163416, 100524, 163479,
                                                                       59304, 59364, 102054,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 165306, 0, 3,
                                                                       163479, 100569, 163542,
                                                                       59364, 59424, 102144,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 165432, 0, 3,
                                                                       163542, 100614, 163605,
                                                                       59424, 59484, 102234,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 165558, 0, 3,
                                                                       163605, 100659, 163668,
                                                                       59484, 59544, 102324,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 165684, 0, 3,
                                                                       163731, 100839, 163794,
                                                                       59664, 59724, 102594,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 165810, 0, 3,
                                                                       163794, 100884, 163857,
                                                                       59724, 59784, 102684,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 165936, 0, 3,
                                                                       163857, 100929, 163920,
                                                                       59784, 59844, 102774,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 166062, 0, 3,
                                                                       163920, 100974, 163983,
                                                                       59844, 59904, 102864,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 166188, 0, 3,
                                                                       163983, 101019, 164046,
                                                                       59904, 59964, 102954,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 166314, 0, 3,
                                                                       164046, 101064, 164109,
                                                                       59964, 60024, 103044,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 166440, 0, 3,
                                                                       164109, 101109, 164172,
                                                                       60024, 60084, 103134,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 166566, 0, 3,
                                                                       164172, 101154, 164235,
                                                                       60084, 60144, 103224,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 166692, 0, 3,
                                                                       164235, 101199, 164298,
                                                                       60144, 60204, 103314,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 166818, 0, 3,
                                                                       164298, 101244, 164361,
                                                                       60204, 60264, 103404,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 166944, 0, 3,
                                                                       164424, 101514, 164550,
                                                                       60384, 60484, 103794,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 167154, 0, 3,
                                                                       164550, 101604, 164676,
                                                                       60484, 60584, 103944,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 167364, 0, 3,
                                                                       164676, 101694, 164802,
                                                                       60584, 60684, 104094,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 167574, 0, 3,
                                                                       164802, 101784, 164928,
                                                                       60684, 60784, 104244,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 167784, 0, 3,
                                                                       164928, 101874, 165054,
                                                                       60784, 60884, 104394,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 167994, 0, 3,
                                                                       165054, 101964, 165180,
                                                                       60884, 60984, 104544,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 168204, 0, 3,
                                                                       165180, 102054, 165306,
                                                                       60984, 61084, 104694,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 168414, 0, 3,
                                                                       165306, 102144, 165432,
                                                                       61084, 61184, 104844,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 168624, 0, 3,
                                                                       165432, 102234, 165558,
                                                                       61184, 61284, 104994,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 168834, 0, 3,
                                                                       165684, 102594, 165810,
                                                                       61484, 61584, 105444,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 169044, 0, 3,
                                                                       165810, 102684, 165936,
                                                                       61584, 61684, 105594,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 169254, 0, 3,
                                                                       165936, 102774, 166062,
                                                                       61684, 61784, 105744,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 169464, 0, 3,
                                                                       166062, 102864, 166188,
                                                                       61784, 61884, 105894,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 169674, 0, 3,
                                                                       166188, 102954, 166314,
                                                                       61884, 61984, 106044,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 169884, 0, 3,
                                                                       166314, 103044, 166440,
                                                                       61984, 62084, 106194,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 170094, 0, 3,
                                                                       166440, 103134, 166566,
                                                                       62084, 62184, 106344,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 170304, 0, 3,
                                                                       166566, 103224, 166692,
                                                                       62184, 62284, 106494,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 170514, 0, 3,
                                                                       166692, 103314, 166818,
                                                                       62284, 62384, 106644,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 170724, 0, 3,
                                                                       166944, 103794, 167154,
                                                                       62584, 62734, 107244,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 171039, 0, 3,
                                                                       167154, 103944, 167364,
                                                                       62734, 62884, 107469,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 171354, 0, 3,
                                                                       167364, 104094, 167574,
                                                                       62884, 63034, 107694,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 171669, 0, 3,
                                                                       167574, 104244, 167784,
                                                                       63034, 63184, 107919,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 171984, 0, 3,
                                                                       167784, 104394, 167994,
                                                                       63184, 63334, 108144,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 172299, 0, 3,
                                                                       167994, 104544, 168204,
                                                                       63334, 63484, 108369,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 172614, 0, 3,
                                                                       168204, 104694, 168414,
                                                                       63484, 63634, 108594,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 172929, 0, 3,
                                                                       168414, 104844, 168624,
                                                                       63634, 63784, 108819,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 173244, 0, 3,
                                                                       168834, 105444, 169044,
                                                                       64084, 64234, 109494,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 173559, 0, 3,
                                                                       169044, 105594, 169254,
                                                                       64234, 64384, 109719,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 173874, 0, 3,
                                                                       169254, 105744, 169464,
                                                                       64384, 64534, 109944,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 174189, 0, 3,
                                                                       169464, 105894, 169674,
                                                                       64534, 64684, 110169,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 174504, 0, 3,
                                                                       169674, 106044, 169884,
                                                                       64684, 64834, 110394,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 174819, 0, 3,
                                                                       169884, 106194, 170094,
                                                                       64834, 64984, 110619,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 175134, 0, 3,
                                                                       170094, 106344, 170304,
                                                                       64984, 65134, 110844,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 175449, 0, 3,
                                                                       170304, 106494, 170514,
                                                                       65134, 65284, 111069,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 175764, 0, 3,
                                                                       170724, 107244, 171039,
                                                                       65584, 65794, 111924,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 176205, 0, 3,
                                                                       171039, 107469, 171354,
                                                                       65794, 66004, 112239,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 176646, 0, 3,
                                                                       171354, 107694, 171669,
                                                                       66004, 66214, 112554,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 177087, 0, 3,
                                                                       171669, 107919, 171984,
                                                                       66214, 66424, 112869,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 177528, 0, 3,
                                                                       171984, 108144, 172299,
                                                                       66424, 66634, 113184,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 177969, 0, 3,
                                                                       172299, 108369, 172614,
                                                                       66634, 66844, 113499,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 178410, 0, 3,
                                                                       172614, 108594, 172929,
                                                                       66844, 67054, 113814,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 178851, 0, 3,
                                                                       173244, 109494, 173559,
                                                                       67474, 67684, 114759,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 179292, 0, 3,
                                                                       173559, 109719, 173874,
                                                                       67684, 67894, 115074,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 179733, 0, 3,
                                                                       173874, 109944, 174189,
                                                                       67894, 68104, 115389,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 180174, 0, 3,
                                                                       174189, 110169, 174504,
                                                                       68104, 68314, 115704,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 180615, 0, 3,
                                                                       174504, 110394, 174819,
                                                                       68314, 68524, 116019,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 181056, 0, 3,
                                                                       174819, 110619, 175134,
                                                                       68524, 68734, 116334,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 181497, 0, 3,
                                                                       175134, 110844, 175449,
                                                                       68734, 68944, 116649,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 181938, 0, 3,
                                                                       175764, 111924, 176205,
                                                                       69364, 69644, 117804,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 182526, 0, 3,
                                                                       176205, 112239, 176646,
                                                                       69644, 69924, 118224,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 183114, 0, 3,
                                                                       176646, 112554, 177087,
                                                                       69924, 70204, 118644,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 183702, 0, 3,
                                                                       177087, 112869, 177528,
                                                                       70204, 70484, 119064,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 184290, 0, 3,
                                                                       177528, 113184, 177969,
                                                                       70484, 70764, 119484,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 184878, 0, 3,
                                                                       177969, 113499, 178410,
                                                                       70764, 71044, 119904,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 185466, 0, 3,
                                                                       178851, 114759, 179292,
                                                                       71604, 71884, 121164,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 186054, 0, 3,
                                                                       179292, 115074, 179733,
                                                                       71884, 72164, 121584,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 186642, 0, 3,
                                                                       179733, 115389, 180174,
                                                                       72164, 72444, 122004,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 187230, 0, 3,
                                                                       180174, 115704, 180615,
                                                                       72444, 72724, 122424,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 187818, 0, 3,
                                                                       180615, 116019, 181056,
                                                                       72724, 73004, 122844,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 188406, 0, 3,
                                                                       181056, 116334, 181497,
                                                                       73004, 73284, 123264,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 188994, 0, 3,
                                                                       181938, 117804, 182526,
                                                                       73844, 74204, 124764,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 189750, 0, 3,
                                                                       182526, 118224, 183114,
                                                                       74204, 74564, 125304,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 190506, 0, 3,
                                                                       183114, 118644, 183702,
                                                                       74564, 74924, 125844,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 191262, 0, 3,
                                                                       183702, 119064, 184290,
                                                                       74924, 75284, 126384,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 192018, 0, 3,
                                                                       184290, 119484, 184878,
                                                                       75284, 75644, 126924,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 192774, 0, 3,
                                                                       185466, 121164, 186054,
                                                                       76364, 76724, 128544,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 193530, 0, 3,
                                                                       186054, 121584, 186642,
                                                                       76724, 77084, 129084,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 194286, 0, 3,
                                                                       186642, 122004, 187230,
                                                                       77084, 77444, 129624,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 195042, 0, 3,
                                                                       187230, 122424, 187818,
                                                                       77444, 77804, 130164,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 195798, 0, 3,
                                                                       187818, 122844, 188406,
                                                                       77804, 78164, 130704,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 196554, 0, 3,
                                                                       188994, 124764, 189750,
                                                                       78884, 79334, 132594,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 197499, 0, 3,
                                                                       189750, 125304, 190506,
                                                                       79334, 79784, 133269,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 198444, 0, 3,
                                                                       190506, 125844, 191262,
                                                                       79784, 80234, 133944,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 199389, 0, 3,
                                                                       191262, 126384, 192018,
                                                                       80234, 80684, 134619,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 200334, 0, 3,
                                                                       192774, 128544, 193530,
                                                                       81584, 82034, 136644,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 201279, 0, 3,
                                                                       193530, 129084, 194286,
                                                                       82034, 82484, 137319,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 202224, 0, 3,
                                                                       194286, 129624, 195042,
                                                                       82484, 82934, 137994,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 203169, 0, 3,
                                                                       195042, 130164, 195798,
                                                                       82934, 83384, 138669,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 204114, 0, 3,
                                                                       196554, 132594, 197499,
                                                                       84284, 84834, 140994,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 205269, 0, 3,
                                                                       197499, 133269, 198444,
                                                                       84834, 85384, 141819,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 206424, 0, 3,
                                                                       198444, 133944, 199389,
                                                                       85384, 85934, 142644,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 207579, 0, 3,
                                                                       200334, 136644, 201279,
                                                                       87034, 87584, 145119,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 208734, 0, 3,
                                                                       201279, 137319, 202224,
                                                                       87584, 88134, 145944,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 209889, 0, 3,
                                                                       202224, 137994, 203169,
                                                                       88134, 88684, 146769,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 211044, 0, 3,
                                                                       204114, 140994, 205269,
                                                                       89784, 90444, 149574,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 212430, 0, 3,
                                                                       205269, 141819, 206424,
                                                                       90444, 91104, 150564,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 213816, 0, 3,
                                                                       207579, 145119, 208734,
                                                                       92424, 93084, 153534,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 215202, 0, 3,
                                                                       208734, 145944, 209889,
                                                                       93084, 93744, 154524,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 216588, 0, 3,
                                                                       211044, 149574, 212430,
                                                                       95064, 95844, 157854,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 218226, 0, 3,
                                                                       213816, 153534, 215202,
                                                                       97404, 98184, 161364,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 219864, 3, 99744,
                                                                       99759, 162534, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 219892, 3, 99759,
                                                                       99774, 162555, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 219920, 3, 99774,
                                                                       99789, 162576, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 219948, 3, 99789,
                                                                       99804, 162597, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 219976, 3, 99804,
                                                                       99819, 162618, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 220004, 3, 99819,
                                                                       99834, 162639, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 220032, 3, 99834,
                                                                       99849, 162660, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 220060, 3, 99849,
                                                                       99864, 162681, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 220088, 3, 99864,
                                                                       99879, 162702, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 220116, 3, 99879,
                                                                       99894, 162723, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 220144, 3, 99894,
                                                                       99909, 162744, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 220172, 3, 99909,
                                                                       99924, 162765, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 220200, 3, 99954,
                                                                       99969, 162786, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 220228, 3, 99969,
                                                                       99984, 162807, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 220256, 3, 99984,
                                                                       99999, 162828, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 220284, 3, 99999,
                                                                       100014, 162849, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 220312, 3, 100014,
                                                                       100029, 162870, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 220340, 3, 100029,
                                                                       100044, 162891, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 220368, 3, 100044,
                                                                       100059, 162912, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 220396, 3, 100059,
                                                                       100074, 162933, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 220424, 3, 100074,
                                                                       100089, 162954, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 220452, 3, 100089,
                                                                       100104, 162975, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 220480, 3, 100104,
                                                                       100119, 162996, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 220508, 3, 100119,
                                                                       100134, 163017, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 220536, 0, 3,
                                                                       219864, 162534, 219892,
                                                                       100164, 100209, 163038,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 220620, 0, 3,
                                                                       219892, 162555, 219920,
                                                                       100209, 100254, 163101,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 220704, 0, 3,
                                                                       219920, 162576, 219948,
                                                                       100254, 100299, 163164,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 220788, 0, 3,
                                                                       219948, 162597, 219976,
                                                                       100299, 100344, 163227,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 220872, 0, 3,
                                                                       219976, 162618, 220004,
                                                                       100344, 100389, 163290,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 220956, 0, 3,
                                                                       220004, 162639, 220032,
                                                                       100389, 100434, 163353,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 221040, 0, 3,
                                                                       220032, 162660, 220060,
                                                                       100434, 100479, 163416,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 221124, 0, 3,
                                                                       220060, 162681, 220088,
                                                                       100479, 100524, 163479,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 221208, 0, 3,
                                                                       220088, 162702, 220116,
                                                                       100524, 100569, 163542,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 221292, 0, 3,
                                                                       220116, 162723, 220144,
                                                                       100569, 100614, 163605,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 221376, 0, 3,
                                                                       220144, 162744, 220172,
                                                                       100614, 100659, 163668,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 221460, 0, 3,
                                                                       220200, 162786, 220228,
                                                                       100749, 100794, 163731,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 221544, 0, 3,
                                                                       220228, 162807, 220256,
                                                                       100794, 100839, 163794,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 221628, 0, 3,
                                                                       220256, 162828, 220284,
                                                                       100839, 100884, 163857,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 221712, 0, 3,
                                                                       220284, 162849, 220312,
                                                                       100884, 100929, 163920,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 221796, 0, 3,
                                                                       220312, 162870, 220340,
                                                                       100929, 100974, 163983,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 221880, 0, 3,
                                                                       220340, 162891, 220368,
                                                                       100974, 101019, 164046,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 221964, 0, 3,
                                                                       220368, 162912, 220396,
                                                                       101019, 101064, 164109,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 222048, 0, 3,
                                                                       220396, 162933, 220424,
                                                                       101064, 101109, 164172,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 222132, 0, 3,
                                                                       220424, 162954, 220452,
                                                                       101109, 101154, 164235,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 222216, 0, 3,
                                                                       220452, 162975, 220480,
                                                                       101154, 101199, 164298,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 222300, 0, 3,
                                                                       220480, 162996, 220508,
                                                                       101199, 101244, 164361,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 222384, 0, 3,
                                                                       220536, 163038, 220620,
                                                                       101334, 101424, 164424,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 222552, 0, 3,
                                                                       220620, 163101, 220704,
                                                                       101424, 101514, 164550,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 222720, 0, 3,
                                                                       220704, 163164, 220788,
                                                                       101514, 101604, 164676,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 222888, 0, 3,
                                                                       220788, 163227, 220872,
                                                                       101604, 101694, 164802,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 223056, 0, 3,
                                                                       220872, 163290, 220956,
                                                                       101694, 101784, 164928,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 223224, 0, 3,
                                                                       220956, 163353, 221040,
                                                                       101784, 101874, 165054,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 223392, 0, 3,
                                                                       221040, 163416, 221124,
                                                                       101874, 101964, 165180,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 223560, 0, 3,
                                                                       221124, 163479, 221208,
                                                                       101964, 102054, 165306,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 223728, 0, 3,
                                                                       221208, 163542, 221292,
                                                                       102054, 102144, 165432,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 223896, 0, 3,
                                                                       221292, 163605, 221376,
                                                                       102144, 102234, 165558,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 224064, 0, 3,
                                                                       221460, 163731, 221544,
                                                                       102414, 102504, 165684,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 224232, 0, 3,
                                                                       221544, 163794, 221628,
                                                                       102504, 102594, 165810,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 224400, 0, 3,
                                                                       221628, 163857, 221712,
                                                                       102594, 102684, 165936,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 224568, 0, 3,
                                                                       221712, 163920, 221796,
                                                                       102684, 102774, 166062,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 224736, 0, 3,
                                                                       221796, 163983, 221880,
                                                                       102774, 102864, 166188,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 224904, 0, 3,
                                                                       221880, 164046, 221964,
                                                                       102864, 102954, 166314,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 225072, 0, 3,
                                                                       221964, 164109, 222048,
                                                                       102954, 103044, 166440,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 225240, 0, 3,
                                                                       222048, 164172, 222132,
                                                                       103044, 103134, 166566,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 225408, 0, 3,
                                                                       222132, 164235, 222216,
                                                                       103134, 103224, 166692,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 225576, 0, 3,
                                                                       222216, 164298, 222300,
                                                                       103224, 103314, 166818,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 225744, 0, 3,
                                                                       222384, 164424, 222552,
                                                                       103494, 103644, 166944,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 226024, 0, 3,
                                                                       222552, 164550, 222720,
                                                                       103644, 103794, 167154,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 226304, 0, 3,
                                                                       222720, 164676, 222888,
                                                                       103794, 103944, 167364,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 226584, 0, 3,
                                                                       222888, 164802, 223056,
                                                                       103944, 104094, 167574,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 226864, 0, 3,
                                                                       223056, 164928, 223224,
                                                                       104094, 104244, 167784,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 227144, 0, 3,
                                                                       223224, 165054, 223392,
                                                                       104244, 104394, 167994,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 227424, 0, 3,
                                                                       223392, 165180, 223560,
                                                                       104394, 104544, 168204,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 227704, 0, 3,
                                                                       223560, 165306, 223728,
                                                                       104544, 104694, 168414,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 227984, 0, 3,
                                                                       223728, 165432, 223896,
                                                                       104694, 104844, 168624,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 228264, 0, 3,
                                                                       224064, 165684, 224232,
                                                                       105144, 105294, 168834,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 228544, 0, 3,
                                                                       224232, 165810, 224400,
                                                                       105294, 105444, 169044,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 228824, 0, 3,
                                                                       224400, 165936, 224568,
                                                                       105444, 105594, 169254,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 229104, 0, 3,
                                                                       224568, 166062, 224736,
                                                                       105594, 105744, 169464,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 229384, 0, 3,
                                                                       224736, 166188, 224904,
                                                                       105744, 105894, 169674,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 229664, 0, 3,
                                                                       224904, 166314, 225072,
                                                                       105894, 106044, 169884,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 229944, 0, 3,
                                                                       225072, 166440, 225240,
                                                                       106044, 106194, 170094,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 230224, 0, 3,
                                                                       225240, 166566, 225408,
                                                                       106194, 106344, 170304,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 230504, 0, 3,
                                                                       225408, 166692, 225576,
                                                                       106344, 106494, 170514,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 230784, 0, 3,
                                                                       225744, 166944, 226024,
                                                                       106794, 107019, 170724,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 231204, 0, 3,
                                                                       226024, 167154, 226304,
                                                                       107019, 107244, 171039,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 231624, 0, 3,
                                                                       226304, 167364, 226584,
                                                                       107244, 107469, 171354,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 232044, 0, 3,
                                                                       226584, 167574, 226864,
                                                                       107469, 107694, 171669,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 232464, 0, 3,
                                                                       226864, 167784, 227144,
                                                                       107694, 107919, 171984,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 232884, 0, 3,
                                                                       227144, 167994, 227424,
                                                                       107919, 108144, 172299,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 233304, 0, 3,
                                                                       227424, 168204, 227704,
                                                                       108144, 108369, 172614,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 233724, 0, 3,
                                                                       227704, 168414, 227984,
                                                                       108369, 108594, 172929,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 234144, 0, 3,
                                                                       228264, 168834, 228544,
                                                                       109044, 109269, 173244,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 234564, 0, 3,
                                                                       228544, 169044, 228824,
                                                                       109269, 109494, 173559,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 234984, 0, 3,
                                                                       228824, 169254, 229104,
                                                                       109494, 109719, 173874,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 235404, 0, 3,
                                                                       229104, 169464, 229384,
                                                                       109719, 109944, 174189,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 235824, 0, 3,
                                                                       229384, 169674, 229664,
                                                                       109944, 110169, 174504,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 236244, 0, 3,
                                                                       229664, 169884, 229944,
                                                                       110169, 110394, 174819,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 236664, 0, 3,
                                                                       229944, 170094, 230224,
                                                                       110394, 110619, 175134,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 237084, 0, 3,
                                                                       230224, 170304, 230504,
                                                                       110619, 110844, 175449,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 237504, 0, 3,
                                                                       230784, 170724, 231204,
                                                                       111294, 111609, 175764,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 238092, 0, 3,
                                                                       231204, 171039, 231624,
                                                                       111609, 111924, 176205,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 238680, 0, 3,
                                                                       231624, 171354, 232044,
                                                                       111924, 112239, 176646,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 239268, 0, 3,
                                                                       232044, 171669, 232464,
                                                                       112239, 112554, 177087,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 239856, 0, 3,
                                                                       232464, 171984, 232884,
                                                                       112554, 112869, 177528,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 240444, 0, 3,
                                                                       232884, 172299, 233304,
                                                                       112869, 113184, 177969,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 241032, 0, 3,
                                                                       233304, 172614, 233724,
                                                                       113184, 113499, 178410,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 241620, 0, 3,
                                                                       234144, 173244, 234564,
                                                                       114129, 114444, 178851,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 242208, 0, 3,
                                                                       234564, 173559, 234984,
                                                                       114444, 114759, 179292,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 242796, 0, 3,
                                                                       234984, 173874, 235404,
                                                                       114759, 115074, 179733,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 243384, 0, 3,
                                                                       235404, 174189, 235824,
                                                                       115074, 115389, 180174,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 243972, 0, 3,
                                                                       235824, 174504, 236244,
                                                                       115389, 115704, 180615,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 244560, 0, 3,
                                                                       236244, 174819, 236664,
                                                                       115704, 116019, 181056,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 245148, 0, 3,
                                                                       236664, 175134, 237084,
                                                                       116019, 116334, 181497,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 245736, 0, 3,
                                                                       237504, 175764, 238092,
                                                                       116964, 117384, 181938,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 246520, 0, 3,
                                                                       238092, 176205, 238680,
                                                                       117384, 117804, 182526,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 247304, 0, 3,
                                                                       238680, 176646, 239268,
                                                                       117804, 118224, 183114,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 248088, 0, 3,
                                                                       239268, 177087, 239856,
                                                                       118224, 118644, 183702,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 248872, 0, 3,
                                                                       239856, 177528, 240444,
                                                                       118644, 119064, 184290,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 249656, 0, 3,
                                                                       240444, 177969, 241032,
                                                                       119064, 119484, 184878,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 250440, 0, 3,
                                                                       241620, 178851, 242208,
                                                                       120324, 120744, 185466,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 251224, 0, 3,
                                                                       242208, 179292, 242796,
                                                                       120744, 121164, 186054,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 252008, 0, 3,
                                                                       242796, 179733, 243384,
                                                                       121164, 121584, 186642,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 252792, 0, 3,
                                                                       243384, 180174, 243972,
                                                                       121584, 122004, 187230,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 253576, 0, 3,
                                                                       243972, 180615, 244560,
                                                                       122004, 122424, 187818,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 254360, 0, 3,
                                                                       244560, 181056, 245148,
                                                                       122424, 122844, 188406,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 255144, 0, 3,
                                                                       245736, 181938, 246520,
                                                                       123684, 124224, 188994,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 256152, 0, 3,
                                                                       246520, 182526, 247304,
                                                                       124224, 124764, 189750,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 257160, 0, 3,
                                                                       247304, 183114, 248088,
                                                                       124764, 125304, 190506,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 258168, 0, 3,
                                                                       248088, 183702, 248872,
                                                                       125304, 125844, 191262,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 259176, 0, 3,
                                                                       248872, 184290, 249656,
                                                                       125844, 126384, 192018,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 260184, 0, 3,
                                                                       250440, 185466, 251224,
                                                                       127464, 128004, 192774,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 261192, 0, 3,
                                                                       251224, 186054, 252008,
                                                                       128004, 128544, 193530,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 262200, 0, 3,
                                                                       252008, 186642, 252792,
                                                                       128544, 129084, 194286,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 263208, 0, 3,
                                                                       252792, 187230, 253576,
                                                                       129084, 129624, 195042,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 264216, 0, 3,
                                                                       253576, 187818, 254360,
                                                                       129624, 130164, 195798,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 265224, 0, 3,
                                                                       255144, 188994, 256152,
                                                                       131244, 131919, 196554,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 266484, 0, 3,
                                                                       256152, 189750, 257160,
                                                                       131919, 132594, 197499,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 267744, 0, 3,
                                                                       257160, 190506, 258168,
                                                                       132594, 133269, 198444,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 269004, 0, 3,
                                                                       258168, 191262, 259176,
                                                                       133269, 133944, 199389,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 270264, 0, 3,
                                                                       260184, 192774, 261192,
                                                                       135294, 135969, 200334,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 271524, 0, 3,
                                                                       261192, 193530, 262200,
                                                                       135969, 136644, 201279,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 272784, 0, 3,
                                                                       262200, 194286, 263208,
                                                                       136644, 137319, 202224,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 274044, 0, 3,
                                                                       263208, 195042, 264216,
                                                                       137319, 137994, 203169,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 275304, 0, 3,
                                                                       265224, 196554, 266484,
                                                                       139344, 140169, 204114,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 276844, 0, 3,
                                                                       266484, 197499, 267744,
                                                                       140169, 140994, 205269,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 278384, 0, 3,
                                                                       267744, 198444, 269004,
                                                                       140994, 141819, 206424,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 279924, 0, 3,
                                                                       270264, 200334, 271524,
                                                                       143469, 144294, 207579,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 281464, 0, 3,
                                                                       271524, 201279, 272784,
                                                                       144294, 145119, 208734,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 283004, 0, 3,
                                                                       272784, 202224, 274044,
                                                                       145119, 145944, 209889,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 284544, 0, 3,
                                                                       275304, 204114, 276844,
                                                                       147594, 148584, 211044,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 286392, 0, 3,
                                                                       276844, 205269, 278384,
                                                                       148584, 149574, 212430,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 288240, 0, 3,
                                                                       279924, 207579, 281464,
                                                                       151554, 152544, 213816,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 290088, 0, 3,
                                                                       281464, 208734, 283004,
                                                                       152544, 153534, 215202,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 291936, 0, 3,
                                                                       284544, 211044, 286392,
                                                                       155514, 156684, 216588,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 294120, 0, 3,
                                                                       288240, 213816, 290088,
                                                                       159024, 160194, 218226,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 296304, 245736, 784, ncols);

                    simdfunc::contract_primitives(buffer, 297452, 250440, 784, ncols);

                    simdfunc::contract_primitives(buffer, 298600, 255144, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 300076, 260184, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 301552, 265224, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 303397, 270264, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 305242, 275304, 1540, ncols);

                    simdfunc::contract_primitives(buffer, 307497, 279924, 1540, ncols);

                    simdfunc::contract_primitives(buffer, 309752, 284544, 1848, ncols);

                    simdfunc::contract_primitives(buffer, 312458, 288240, 1848, ncols);

                    simdfunc::contract_primitives(buffer, 315164, 291936, 2184, ncols);

                    simdfunc::contract_primitives(buffer, 318362, 294120, 2184, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 297088, 296304, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 298236, 297452, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 299608, 298600, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 301084, 300076, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 302812, 301552, 45, 1, nmax);

        simdtrf::transform_i_inner(buffer, 304657, 303397, 45, 1, nmax);

        simdtrf::transform_i_inner(buffer, 306782, 305242, 55, 1, nmax);

        simdtrf::transform_i_inner(buffer, 309037, 307497, 55, 1, nmax);

        simdtrf::transform_i_inner(buffer, 311600, 309752, 66, 1, nmax);

        simdtrf::transform_i_inner(buffer, 314306, 312458, 66, 1, nmax);

        simdtrf::transform_i_inner(buffer, 317348, 315164, 78, 1, nmax);

        simdtrf::transform_i_inner(buffer, 320546, 318362, 78, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 321560, 297088, 299608, 13,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 322652, 298236, 301084, 13,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 323744, 299608, 302812, 13,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 325148, 301084, 304657, 13,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 326552, 302812, 306782, 13,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 328307, 304657, 309037, 13,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 330062, 306782, 311600, 13,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 332207, 309037, 314306, 13,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 334352, 311600, 317348, 13,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 336926, 314306, 320546, 13,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 339500, 321560, 323744, 13,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 341684, 322652, 325148, 13,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 343868, 323744, 326552, 13,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 346676, 325148, 328307, 13,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 349484, 326552, 330062, 13,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 352994, 328307, 332207, 13,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 356504, 330062, 334352, 13,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 360794, 332207, 336926, 13,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 365084, 339500, 343868, 13,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 368724, 341684, 346676, 13,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 372364, 343868, 349484, 13,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 377044, 346676, 352994, 13,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 381724, 349484, 356504, 13,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 387574, 352994, 360794, 13,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 393424, 365084, 372364, 13,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 398884, 368724, 377044, 13,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 404344, 372364, 381724, 13,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 411364, 377044, 387574, 13,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 418384, 393424, 404344, 13,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 426028, 398884, 411364, 13,
                                             nmax);

        simdtrf::transform_h_inner(buffer, 433672, 426028, 28, 13, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 433672, 143, nmax);

        simdtrf::transform_h_inner(buffer, 433672, 418384, 28, 13, nmax);

        simdtrf::transform_i_outer(values + 1859 * nvalues + n * npairs, nvalues, buffer, 433672,
                                   143, nmax);
    }

    for (size_t m = 0; m < 3718; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
