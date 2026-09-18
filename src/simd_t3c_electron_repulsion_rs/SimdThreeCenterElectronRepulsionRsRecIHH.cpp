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


#include "SimdThreeCenterElectronRepulsionRsRecIHH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
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
compute_rs_ihh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_ihh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 305796, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 3146 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 305796, 187832, 18854, dimensions);

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

                    simdfunc::compute_t3c_erf_boys_function(buffer, coordinates, 6, 3, {1, 2, 3,
                                                            4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                                                            15, 16}, ncols, fj, i * nprim_b + j,
                                                            fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 23, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                                        16}, ncols, fj, i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 40, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 43, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 46, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 49, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 52, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 55, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 58, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 61, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 64, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 67, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 70, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 73, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 76, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 79, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 82, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 85, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 88, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 91, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 94, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 97, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 100, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 103, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 106, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 109, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 112, 0, 3, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 115, 0, 3, 34, 35,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 118, 0, 3, 35, 36,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 121, 0, 3, 36, 37,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 124, 0, 3, 37, 38,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 127, 0, 3, 38, 39,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 130, 0, 3, 7, 8,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 136, 0, 3, 8, 9,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 142, 0, 3, 9, 10,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 148, 0, 3, 10, 11,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 154, 0, 3, 11, 12,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 160, 0, 3, 12, 13,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 166, 0, 3, 13, 14,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 172, 0, 3, 14, 15,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 178, 0, 3, 15, 16,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 184, 0, 3, 16, 17,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 190, 0, 3, 17, 18,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 196, 0, 3, 18, 19,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 202, 0, 3, 19, 20,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 208, 0, 3, 20, 21,
                                                                       79, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 214, 0, 3, 24, 25,
                                                                       85, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 220, 0, 3, 25, 26,
                                                                       88, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 226, 0, 3, 26, 27,
                                                                       91, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 232, 0, 3, 27, 28,
                                                                       94, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 238, 0, 3, 28, 29,
                                                                       97, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 244, 0, 3, 29, 30,
                                                                       100, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 250, 0, 3, 30, 31,
                                                                       103, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 256, 0, 3, 31, 32,
                                                                       106, 109, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 262, 0, 3, 32, 33,
                                                                       109, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 268, 0, 3, 33, 34,
                                                                       112, 115, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 274, 0, 3, 34, 35,
                                                                       115, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 280, 0, 3, 35, 36,
                                                                       118, 121, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 286, 0, 3, 36, 37,
                                                                       121, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 292, 0, 3, 37, 38,
                                                                       124, 127, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 298, 0, 3, 40, 43,
                                                                       130, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 308, 0, 3, 43, 46,
                                                                       136, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 318, 0, 3, 46, 49,
                                                                       142, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 328, 0, 3, 49, 52,
                                                                       148, 154, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 338, 0, 3, 52, 55,
                                                                       154, 160, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 348, 0, 3, 55, 58,
                                                                       160, 166, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 58, 61,
                                                                       166, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 368, 0, 3, 61, 64,
                                                                       172, 178, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 378, 0, 3, 64, 67,
                                                                       178, 184, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 388, 0, 3, 67, 70,
                                                                       184, 190, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 398, 0, 3, 70, 73,
                                                                       190, 196, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 408, 0, 3, 73, 76,
                                                                       196, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 418, 0, 3, 76, 79,
                                                                       202, 208, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 85, 88,
                                                                       214, 220, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 438, 0, 3, 88, 91,
                                                                       220, 226, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 448, 0, 3, 91, 94,
                                                                       226, 232, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 458, 0, 3, 94, 97,
                                                                       232, 238, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 468, 0, 3, 97,
                                                                       100, 238, 244, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 478, 0, 3, 100,
                                                                       103, 244, 250, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 488, 0, 3, 103,
                                                                       106, 250, 256, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 498, 0, 3, 106,
                                                                       109, 256, 262, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 508, 0, 3, 109,
                                                                       112, 262, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 518, 0, 3, 112,
                                                                       115, 268, 274, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 528, 0, 3, 115,
                                                                       118, 274, 280, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 538, 0, 3, 118,
                                                                       121, 280, 286, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 548, 0, 3, 121,
                                                                       124, 286, 292, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 558, 0, 3, 130,
                                                                       136, 298, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 573, 0, 3, 136,
                                                                       142, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 588, 0, 3, 142,
                                                                       148, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 603, 0, 3, 148,
                                                                       154, 328, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 618, 0, 3, 154,
                                                                       160, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 633, 0, 3, 160,
                                                                       166, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 648, 0, 3, 166,
                                                                       172, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 663, 0, 3, 172,
                                                                       178, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 678, 0, 3, 178,
                                                                       184, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 693, 0, 3, 184,
                                                                       190, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 708, 0, 3, 190,
                                                                       196, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 723, 0, 3, 196,
                                                                       202, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 738, 0, 3, 214,
                                                                       220, 428, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 753, 0, 3, 220,
                                                                       226, 438, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 768, 0, 3, 226,
                                                                       232, 448, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 783, 0, 3, 232,
                                                                       238, 458, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 798, 0, 3, 238,
                                                                       244, 468, 478, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 813, 0, 3, 244,
                                                                       250, 478, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 828, 0, 3, 250,
                                                                       256, 488, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 843, 0, 3, 256,
                                                                       262, 498, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 858, 0, 3, 262,
                                                                       268, 508, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 873, 0, 3, 268,
                                                                       274, 518, 528, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 888, 0, 3, 274,
                                                                       280, 528, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 903, 0, 3, 280,
                                                                       286, 538, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 918, 0, 3, 298,
                                                                       308, 558, 573, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 939, 0, 3, 308,
                                                                       318, 573, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 960, 0, 3, 318,
                                                                       328, 588, 603, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 981, 0, 3, 328,
                                                                       338, 603, 618, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1002, 0, 3, 338,
                                                                       348, 618, 633, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1023, 0, 3, 348,
                                                                       358, 633, 648, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1044, 0, 3, 358,
                                                                       368, 648, 663, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1065, 0, 3, 368,
                                                                       378, 663, 678, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1086, 0, 3, 378,
                                                                       388, 678, 693, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1107, 0, 3, 388,
                                                                       398, 693, 708, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1128, 0, 3, 398,
                                                                       408, 708, 723, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1149, 0, 3, 428,
                                                                       438, 738, 753, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1170, 0, 3, 438,
                                                                       448, 753, 768, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1191, 0, 3, 448,
                                                                       458, 768, 783, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1212, 0, 3, 458,
                                                                       468, 783, 798, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1233, 0, 3, 468,
                                                                       478, 798, 813, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1254, 0, 3, 478,
                                                                       488, 813, 828, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1275, 0, 3, 488,
                                                                       498, 828, 843, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1296, 0, 3, 498,
                                                                       508, 843, 858, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1317, 0, 3, 508,
                                                                       518, 858, 873, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1338, 0, 3, 518,
                                                                       528, 873, 888, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1359, 0, 3, 528,
                                                                       538, 888, 903, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 558,
                                                                       573, 918, 939, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1408, 0, 3, 573,
                                                                       588, 939, 960, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1436, 0, 3, 588,
                                                                       603, 960, 981, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1464, 0, 3, 603,
                                                                       618, 981, 1002, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 618,
                                                                       633, 1002, 1023, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 633,
                                                                       648, 1023, 1044, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1548, 0, 3, 648,
                                                                       663, 1044, 1065, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1576, 0, 3, 663,
                                                                       678, 1065, 1086, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1604, 0, 3, 678,
                                                                       693, 1086, 1107, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1632, 0, 3, 693,
                                                                       708, 1107, 1128, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1660, 0, 3, 738,
                                                                       753, 1149, 1170, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 753,
                                                                       768, 1170, 1191, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1716, 0, 3, 768,
                                                                       783, 1191, 1212, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1744, 0, 3, 783,
                                                                       798, 1212, 1233, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1772, 0, 3, 798,
                                                                       813, 1233, 1254, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1800, 0, 3, 813,
                                                                       828, 1254, 1275, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1828, 0, 3, 828,
                                                                       843, 1275, 1296, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1856, 0, 3, 843,
                                                                       858, 1296, 1317, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1884, 0, 3, 858,
                                                                       873, 1317, 1338, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1912, 0, 3, 873,
                                                                       888, 1338, 1359, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1940, 0, 3, 918,
                                                                       939, 1380, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1976, 0, 3, 939,
                                                                       960, 1408, 1436, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2012, 0, 3, 960,
                                                                       981, 1436, 1464, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2048, 0, 3, 981,
                                                                       1002, 1464, 1492, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2084, 0, 3, 1002,
                                                                       1023, 1492, 1520, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2120, 0, 3, 1023,
                                                                       1044, 1520, 1548, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2156, 0, 3, 1044,
                                                                       1065, 1548, 1576, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2192, 0, 3, 1065,
                                                                       1086, 1576, 1604, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2228, 0, 3, 1086,
                                                                       1107, 1604, 1632, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2264, 0, 3, 1149,
                                                                       1170, 1660, 1688, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2300, 0, 3, 1170,
                                                                       1191, 1688, 1716, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2336, 0, 3, 1191,
                                                                       1212, 1716, 1744, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2372, 0, 3, 1212,
                                                                       1233, 1744, 1772, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2408, 0, 3, 1233,
                                                                       1254, 1772, 1800, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2444, 0, 3, 1254,
                                                                       1275, 1800, 1828, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2480, 0, 3, 1275,
                                                                       1296, 1828, 1856, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2516, 0, 3, 1296,
                                                                       1317, 1856, 1884, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2552, 0, 3, 1317,
                                                                       1338, 1884, 1912, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2588, 0, 3, 1380,
                                                                       1408, 1940, 1976, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2633, 0, 3, 1408,
                                                                       1436, 1976, 2012, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2678, 0, 3, 1436,
                                                                       1464, 2012, 2048, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2723, 0, 3, 1464,
                                                                       1492, 2048, 2084, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2768, 0, 3, 1492,
                                                                       1520, 2084, 2120, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2813, 0, 3, 1520,
                                                                       1548, 2120, 2156, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2858, 0, 3, 1548,
                                                                       1576, 2156, 2192, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2903, 0, 3, 1576,
                                                                       1604, 2192, 2228, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2948, 0, 3, 1660,
                                                                       1688, 2264, 2300, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2993, 0, 3, 1688,
                                                                       1716, 2300, 2336, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3038, 0, 3, 1716,
                                                                       1744, 2336, 2372, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3083, 0, 3, 1744,
                                                                       1772, 2372, 2408, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3128, 0, 3, 1772,
                                                                       1800, 2408, 2444, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3173, 0, 3, 1800,
                                                                       1828, 2444, 2480, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3218, 0, 3, 1828,
                                                                       1856, 2480, 2516, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3263, 0, 3, 1856,
                                                                       1884, 2516, 2552, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3308, 0, 3, 1940,
                                                                       1976, 2588, 2633, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3363, 0, 3, 1976,
                                                                       2012, 2633, 2678, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3418, 0, 3, 2012,
                                                                       2048, 2678, 2723, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3473, 0, 3, 2048,
                                                                       2084, 2723, 2768, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3528, 0, 3, 2084,
                                                                       2120, 2768, 2813, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3583, 0, 3, 2120,
                                                                       2156, 2813, 2858, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3638, 0, 3, 2156,
                                                                       2192, 2858, 2903, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3693, 0, 3, 2264,
                                                                       2300, 2948, 2993, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3748, 0, 3, 2300,
                                                                       2336, 2993, 3038, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3803, 0, 3, 2336,
                                                                       2372, 3038, 3083, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3858, 0, 3, 2372,
                                                                       2408, 3083, 3128, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3913, 0, 3, 2408,
                                                                       2444, 3128, 3173, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3968, 0, 3, 2444,
                                                                       2480, 3173, 3218, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4023, 0, 3, 2480,
                                                                       2516, 3218, 3263, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4078, 0, 3, 2588,
                                                                       2633, 3308, 3363, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4144, 0, 3, 2633,
                                                                       2678, 3363, 3418, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4210, 0, 3, 2678,
                                                                       2723, 3418, 3473, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4276, 0, 3, 2723,
                                                                       2768, 3473, 3528, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4342, 0, 3, 2768,
                                                                       2813, 3528, 3583, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4408, 0, 3, 2813,
                                                                       2858, 3583, 3638, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4474, 0, 3, 2948,
                                                                       2993, 3693, 3748, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4540, 0, 3, 2993,
                                                                       3038, 3748, 3803, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4606, 0, 3, 3038,
                                                                       3083, 3803, 3858, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4672, 0, 3, 3083,
                                                                       3128, 3858, 3913, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4738, 0, 3, 3128,
                                                                       3173, 3913, 3968, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4804, 0, 3, 3173,
                                                                       3218, 3968, 4023, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4870, 0, 3, 3308,
                                                                       3363, 4078, 4144, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4948, 0, 3, 3363,
                                                                       3418, 4144, 4210, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 5026, 0, 3, 3418,
                                                                       3473, 4210, 4276, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 5104, 0, 3, 3473,
                                                                       3528, 4276, 4342, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 5182, 0, 3, 3528,
                                                                       3583, 4342, 4408, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 5260, 0, 3, 3693,
                                                                       3748, 4474, 4540, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 5338, 0, 3, 3748,
                                                                       3803, 4540, 4606, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 5416, 0, 3, 3803,
                                                                       3858, 4606, 4672, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 5494, 0, 3, 3858,
                                                                       3913, 4672, 4738, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 5572, 0, 3, 3913,
                                                                       3968, 4738, 4804, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5650, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5653, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5656, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5659, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5662, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5665, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5668, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5671, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5674, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5677, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5680, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5683, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5686, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5689, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5692, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5695, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5698, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5701, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5704, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5707, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5710, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5713, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5716, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5719, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5722, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5725, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5728, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5731, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5734, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5737, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5740, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5743, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5746, 3, 7, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5755, 3, 8, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5764, 3, 9, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5773, 3, 10, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5782, 3, 11, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5791, 3, 12, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5800, 3, 13, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5809, 3, 14, 61,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5818, 3, 15, 64,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5827, 3, 16, 67,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5836, 3, 17, 70,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5845, 3, 18, 73,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5854, 3, 19, 76,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5863, 3, 20, 79,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5872, 3, 21, 82,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5881, 3, 24, 85,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5890, 3, 25, 88,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5899, 3, 26, 91,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5908, 3, 27, 94,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5917, 3, 28, 97,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5926, 3, 29, 100,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5935, 3, 30, 103,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5944, 3, 31, 106,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5953, 3, 32, 109,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5962, 3, 33, 112,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5971, 3, 34, 115,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5980, 3, 35, 118,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5989, 3, 36, 121,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5998, 3, 37, 124,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 6007, 3, 38, 127,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6016, 3, 40, 130,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6034, 3, 43, 136,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6052, 3, 46, 142,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6070, 3, 49, 148,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6088, 3, 52, 154,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6106, 3, 55, 160,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6124, 3, 58, 166,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6142, 3, 61, 172,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6160, 3, 64, 178,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6178, 3, 67, 184,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6196, 3, 70, 190,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6214, 3, 73, 196,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6232, 3, 76, 202,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6250, 3, 79, 208,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6268, 3, 85, 214,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6286, 3, 88, 220,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6304, 3, 91, 226,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6322, 3, 94, 232,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6340, 3, 97, 238,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6358, 3, 100, 244,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6376, 3, 103, 250,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6394, 3, 106, 256,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6412, 3, 109, 262,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6430, 3, 112, 268,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6448, 3, 115, 274,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6466, 3, 118, 280,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6484, 3, 121, 286,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6502, 3, 124, 292,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6520, 3, 130, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6550, 3, 136, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6580, 3, 142, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6610, 3, 148, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6640, 3, 154, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6670, 3, 160, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6700, 3, 166, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6730, 3, 172, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6760, 3, 178, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6790, 3, 184, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6820, 3, 190, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6850, 3, 196, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6880, 3, 202, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6910, 3, 214, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6940, 3, 220, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6970, 3, 226, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 7000, 3, 232, 458,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 7030, 3, 238, 468,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 7060, 3, 244, 478,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 7090, 3, 250, 488,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 7120, 3, 256, 498,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 7150, 3, 262, 508,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 7180, 3, 268, 518,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 7210, 3, 274, 528,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 7240, 3, 280, 538,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 7270, 3, 286, 548,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7300, 3, 298, 558,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7345, 3, 308, 573,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7390, 3, 318, 588,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7435, 3, 328, 603,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7480, 3, 338, 618,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7525, 3, 348, 633,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7570, 3, 358, 648,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7615, 3, 368, 663,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7660, 3, 378, 678,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7705, 3, 388, 693,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7750, 3, 398, 708,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7795, 3, 408, 723,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7840, 3, 428, 738,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7885, 3, 438, 753,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7930, 3, 448, 768,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7975, 3, 458, 783,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8020, 3, 468, 798,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8065, 3, 478, 813,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8110, 3, 488, 828,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8155, 3, 498, 843,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8200, 3, 508, 858,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8245, 3, 518, 873,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8290, 3, 528, 888,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8335, 3, 538, 903,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8380, 3, 558, 918,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8443, 3, 573, 939,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8506, 3, 588, 960,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8569, 3, 603, 981,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8632, 3, 618,
                                                                       1002, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8695, 3, 633,
                                                                       1023, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8758, 3, 648,
                                                                       1044, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8821, 3, 663,
                                                                       1065, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8884, 3, 678,
                                                                       1086, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8947, 3, 693,
                                                                       1107, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9010, 3, 708,
                                                                       1128, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9073, 3, 738,
                                                                       1149, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9136, 3, 753,
                                                                       1170, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9199, 3, 768,
                                                                       1191, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9262, 3, 783,
                                                                       1212, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9325, 3, 798,
                                                                       1233, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9388, 3, 813,
                                                                       1254, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9451, 3, 828,
                                                                       1275, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9514, 3, 843,
                                                                       1296, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9577, 3, 858,
                                                                       1317, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9640, 3, 873,
                                                                       1338, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9703, 3, 888,
                                                                       1359, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9766, 3, 918,
                                                                       1380, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9850, 3, 939,
                                                                       1408, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9934, 3, 960,
                                                                       1436, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10018, 3, 981,
                                                                       1464, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10102, 3, 1002,
                                                                       1492, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10186, 3, 1023,
                                                                       1520, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10270, 3, 1044,
                                                                       1548, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10354, 3, 1065,
                                                                       1576, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10438, 3, 1086,
                                                                       1604, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10522, 3, 1107,
                                                                       1632, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10606, 3, 1149,
                                                                       1660, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10690, 3, 1170,
                                                                       1688, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10774, 3, 1191,
                                                                       1716, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10858, 3, 1212,
                                                                       1744, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10942, 3, 1233,
                                                                       1772, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11026, 3, 1254,
                                                                       1800, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11110, 3, 1275,
                                                                       1828, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11194, 3, 1296,
                                                                       1856, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11278, 3, 1317,
                                                                       1884, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11362, 3, 1338,
                                                                       1912, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11446, 3, 1380,
                                                                       1940, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11554, 3, 1408,
                                                                       1976, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11662, 3, 1436,
                                                                       2012, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11770, 3, 1464,
                                                                       2048, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11878, 3, 1492,
                                                                       2084, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11986, 3, 1520,
                                                                       2120, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12094, 3, 1548,
                                                                       2156, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12202, 3, 1576,
                                                                       2192, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12310, 3, 1604,
                                                                       2228, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12418, 3, 1660,
                                                                       2264, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12526, 3, 1688,
                                                                       2300, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12634, 3, 1716,
                                                                       2336, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12742, 3, 1744,
                                                                       2372, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12850, 3, 1772,
                                                                       2408, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12958, 3, 1800,
                                                                       2444, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13066, 3, 1828,
                                                                       2480, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13174, 3, 1856,
                                                                       2516, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13282, 3, 1884,
                                                                       2552, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13390, 3, 1940,
                                                                       2588, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13525, 3, 1976,
                                                                       2633, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13660, 3, 2012,
                                                                       2678, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13795, 3, 2048,
                                                                       2723, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13930, 3, 2084,
                                                                       2768, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14065, 3, 2120,
                                                                       2813, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14200, 3, 2156,
                                                                       2858, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14335, 3, 2192,
                                                                       2903, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14470, 3, 2264,
                                                                       2948, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14605, 3, 2300,
                                                                       2993, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14740, 3, 2336,
                                                                       3038, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14875, 3, 2372,
                                                                       3083, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15010, 3, 2408,
                                                                       3128, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15145, 3, 2444,
                                                                       3173, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15280, 3, 2480,
                                                                       3218, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 15415, 3, 2516,
                                                                       3263, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15550, 3, 2588,
                                                                       3308, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15715, 3, 2633,
                                                                       3363, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15880, 3, 2678,
                                                                       3418, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16045, 3, 2723,
                                                                       3473, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16210, 3, 2768,
                                                                       3528, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16375, 3, 2813,
                                                                       3583, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16540, 3, 2858,
                                                                       3638, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16705, 3, 2948,
                                                                       3693, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16870, 3, 2993,
                                                                       3748, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17035, 3, 3038,
                                                                       3803, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17200, 3, 3083,
                                                                       3858, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17365, 3, 3128,
                                                                       3913, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17530, 3, 3173,
                                                                       3968, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 17695, 3, 3218,
                                                                       4023, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 17860, 3, 3308,
                                                                       4078, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 18058, 3, 3363,
                                                                       4144, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 18256, 3, 3418,
                                                                       4210, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 18454, 3, 3473,
                                                                       4276, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 18652, 3, 3528,
                                                                       4342, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 18850, 3, 3583,
                                                                       4408, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 19048, 3, 3693,
                                                                       4474, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 19246, 3, 3748,
                                                                       4540, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 19444, 3, 3803,
                                                                       4606, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 19642, 3, 3858,
                                                                       4672, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 19840, 3, 3913,
                                                                       4738, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 20038, 3, 3968,
                                                                       4804, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 20236, 3, 4078,
                                                                       4870, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 20470, 3, 4144,
                                                                       4948, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 20704, 3, 4210,
                                                                       5026, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 20938, 3, 4276,
                                                                       5104, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 21172, 3, 4342,
                                                                       5182, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 21406, 3, 4474,
                                                                       5260, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 21640, 3, 4540,
                                                                       5338, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 21874, 3, 4606,
                                                                       5416, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 22108, 3, 4672,
                                                                       5494, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 22342, 3, 4738,
                                                                       5572, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22576, 3, 7, 8,
                                                                       5656, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22582, 3, 8, 9,
                                                                       5659, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22588, 3, 9, 10,
                                                                       5662, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22594, 3, 10, 11,
                                                                       5665, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22600, 3, 11, 12,
                                                                       5668, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22606, 3, 12, 13,
                                                                       5671, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22612, 3, 13, 14,
                                                                       5674, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22618, 3, 14, 15,
                                                                       5677, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22624, 3, 15, 16,
                                                                       5680, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22630, 3, 16, 17,
                                                                       5683, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22636, 3, 17, 18,
                                                                       5686, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22642, 3, 18, 19,
                                                                       5689, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22648, 3, 19, 20,
                                                                       5692, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22654, 3, 20, 21,
                                                                       5695, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22660, 3, 24, 25,
                                                                       5704, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22666, 3, 25, 26,
                                                                       5707, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22672, 3, 26, 27,
                                                                       5710, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22678, 3, 27, 28,
                                                                       5713, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22684, 3, 28, 29,
                                                                       5716, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22690, 3, 29, 30,
                                                                       5719, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22696, 3, 30, 31,
                                                                       5722, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22702, 3, 31, 32,
                                                                       5725, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22708, 3, 32, 33,
                                                                       5728, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22714, 3, 33, 34,
                                                                       5731, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22720, 3, 34, 35,
                                                                       5734, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22726, 3, 35, 36,
                                                                       5737, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22732, 3, 36, 37,
                                                                       5740, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22738, 3, 37, 38,
                                                                       5743, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22744, 0, 3,
                                                                       22576, 5656, 22582, 5764,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22762, 0, 3,
                                                                       22582, 5659, 22588, 5773,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22780, 0, 3,
                                                                       22588, 5662, 22594, 5782,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22798, 0, 3,
                                                                       22594, 5665, 22600, 5791,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22816, 0, 3,
                                                                       22600, 5668, 22606, 5800,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22834, 0, 3,
                                                                       22606, 5671, 22612, 5809,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22852, 0, 3,
                                                                       22612, 5674, 22618, 5818,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22870, 0, 3,
                                                                       22618, 5677, 22624, 5827,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22888, 0, 3,
                                                                       22624, 5680, 22630, 5836,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22906, 0, 3,
                                                                       22630, 5683, 22636, 5845,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22924, 0, 3,
                                                                       22636, 5686, 22642, 5854,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22942, 0, 3,
                                                                       22642, 5689, 22648, 5863,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22960, 0, 3,
                                                                       22648, 5692, 22654, 5872,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22978, 0, 3,
                                                                       22660, 5704, 22666, 5899,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22996, 0, 3,
                                                                       22666, 5707, 22672, 5908,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 23014, 0, 3,
                                                                       22672, 5710, 22678, 5917,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 23032, 0, 3,
                                                                       22678, 5713, 22684, 5926,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 23050, 0, 3,
                                                                       22684, 5716, 22690, 5935,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 23068, 0, 3,
                                                                       22690, 5719, 22696, 5944,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 23086, 0, 3,
                                                                       22696, 5722, 22702, 5953,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 23104, 0, 3,
                                                                       22702, 5725, 22708, 5962,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 23122, 0, 3,
                                                                       22708, 5728, 22714, 5971,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 23140, 0, 3,
                                                                       22714, 5731, 22720, 5980,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 23158, 0, 3,
                                                                       22720, 5734, 22726, 5989,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 23176, 0, 3,
                                                                       22726, 5737, 22732, 5998,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 23194, 0, 3,
                                                                       22732, 5740, 22738, 6007,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23212, 0, 3,
                                                                       22744, 5764, 22762, 130,
                                                                       136, 6052, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23248, 0, 3,
                                                                       22762, 5773, 22780, 136,
                                                                       142, 6070, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23284, 0, 3,
                                                                       22780, 5782, 22798, 142,
                                                                       148, 6088, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23320, 0, 3,
                                                                       22798, 5791, 22816, 148,
                                                                       154, 6106, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23356, 0, 3,
                                                                       22816, 5800, 22834, 154,
                                                                       160, 6124, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23392, 0, 3,
                                                                       22834, 5809, 22852, 160,
                                                                       166, 6142, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23428, 0, 3,
                                                                       22852, 5818, 22870, 166,
                                                                       172, 6160, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23464, 0, 3,
                                                                       22870, 5827, 22888, 172,
                                                                       178, 6178, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23500, 0, 3,
                                                                       22888, 5836, 22906, 178,
                                                                       184, 6196, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23536, 0, 3,
                                                                       22906, 5845, 22924, 184,
                                                                       190, 6214, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23572, 0, 3,
                                                                       22924, 5854, 22942, 190,
                                                                       196, 6232, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23608, 0, 3,
                                                                       22942, 5863, 22960, 196,
                                                                       202, 6250, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23644, 0, 3,
                                                                       22978, 5899, 22996, 214,
                                                                       220, 6304, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23680, 0, 3,
                                                                       22996, 5908, 23014, 220,
                                                                       226, 6322, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23716, 0, 3,
                                                                       23014, 5917, 23032, 226,
                                                                       232, 6340, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23752, 0, 3,
                                                                       23032, 5926, 23050, 232,
                                                                       238, 6358, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23788, 0, 3,
                                                                       23050, 5935, 23068, 238,
                                                                       244, 6376, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23824, 0, 3,
                                                                       23068, 5944, 23086, 244,
                                                                       250, 6394, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23860, 0, 3,
                                                                       23086, 5953, 23104, 250,
                                                                       256, 6412, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23896, 0, 3,
                                                                       23104, 5962, 23122, 256,
                                                                       262, 6430, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23932, 0, 3,
                                                                       23122, 5971, 23140, 262,
                                                                       268, 6448, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23968, 0, 3,
                                                                       23140, 5980, 23158, 268,
                                                                       274, 6466, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24004, 0, 3,
                                                                       23158, 5989, 23176, 274,
                                                                       280, 6484, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 24040, 0, 3,
                                                                       23176, 5998, 23194, 280,
                                                                       286, 6502, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24076, 0, 3,
                                                                       23212, 6052, 23248, 298,
                                                                       308, 6580, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24136, 0, 3,
                                                                       23248, 6070, 23284, 308,
                                                                       318, 6610, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24196, 0, 3,
                                                                       23284, 6088, 23320, 318,
                                                                       328, 6640, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24256, 0, 3,
                                                                       23320, 6106, 23356, 328,
                                                                       338, 6670, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24316, 0, 3,
                                                                       23356, 6124, 23392, 338,
                                                                       348, 6700, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24376, 0, 3,
                                                                       23392, 6142, 23428, 348,
                                                                       358, 6730, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24436, 0, 3,
                                                                       23428, 6160, 23464, 358,
                                                                       368, 6760, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24496, 0, 3,
                                                                       23464, 6178, 23500, 368,
                                                                       378, 6790, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24556, 0, 3,
                                                                       23500, 6196, 23536, 378,
                                                                       388, 6820, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24616, 0, 3,
                                                                       23536, 6214, 23572, 388,
                                                                       398, 6850, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24676, 0, 3,
                                                                       23572, 6232, 23608, 398,
                                                                       408, 6880, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24736, 0, 3,
                                                                       23644, 6304, 23680, 428,
                                                                       438, 6970, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24796, 0, 3,
                                                                       23680, 6322, 23716, 438,
                                                                       448, 7000, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24856, 0, 3,
                                                                       23716, 6340, 23752, 448,
                                                                       458, 7030, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24916, 0, 3,
                                                                       23752, 6358, 23788, 458,
                                                                       468, 7060, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24976, 0, 3,
                                                                       23788, 6376, 23824, 468,
                                                                       478, 7090, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25036, 0, 3,
                                                                       23824, 6394, 23860, 478,
                                                                       488, 7120, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25096, 0, 3,
                                                                       23860, 6412, 23896, 488,
                                                                       498, 7150, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25156, 0, 3,
                                                                       23896, 6430, 23932, 498,
                                                                       508, 7180, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25216, 0, 3,
                                                                       23932, 6448, 23968, 508,
                                                                       518, 7210, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25276, 0, 3,
                                                                       23968, 6466, 24004, 518,
                                                                       528, 7240, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 25336, 0, 3,
                                                                       24004, 6484, 24040, 528,
                                                                       538, 7270, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25396, 0, 3,
                                                                       24076, 6580, 24136, 558,
                                                                       573, 7390, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25486, 0, 3,
                                                                       24136, 6610, 24196, 573,
                                                                       588, 7435, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25576, 0, 3,
                                                                       24196, 6640, 24256, 588,
                                                                       603, 7480, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25666, 0, 3,
                                                                       24256, 6670, 24316, 603,
                                                                       618, 7525, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25756, 0, 3,
                                                                       24316, 6700, 24376, 618,
                                                                       633, 7570, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25846, 0, 3,
                                                                       24376, 6730, 24436, 633,
                                                                       648, 7615, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25936, 0, 3,
                                                                       24436, 6760, 24496, 648,
                                                                       663, 7660, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26026, 0, 3,
                                                                       24496, 6790, 24556, 663,
                                                                       678, 7705, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26116, 0, 3,
                                                                       24556, 6820, 24616, 678,
                                                                       693, 7750, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26206, 0, 3,
                                                                       24616, 6850, 24676, 693,
                                                                       708, 7795, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26296, 0, 3,
                                                                       24736, 6970, 24796, 738,
                                                                       753, 7930, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26386, 0, 3,
                                                                       24796, 7000, 24856, 753,
                                                                       768, 7975, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26476, 0, 3,
                                                                       24856, 7030, 24916, 768,
                                                                       783, 8020, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26566, 0, 3,
                                                                       24916, 7060, 24976, 783,
                                                                       798, 8065, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26656, 0, 3,
                                                                       24976, 7090, 25036, 798,
                                                                       813, 8110, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26746, 0, 3,
                                                                       25036, 7120, 25096, 813,
                                                                       828, 8155, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26836, 0, 3,
                                                                       25096, 7150, 25156, 828,
                                                                       843, 8200, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 26926, 0, 3,
                                                                       25156, 7180, 25216, 843,
                                                                       858, 8245, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27016, 0, 3,
                                                                       25216, 7210, 25276, 858,
                                                                       873, 8290, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 27106, 0, 3,
                                                                       25276, 7240, 25336, 873,
                                                                       888, 8335, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27196, 0, 3,
                                                                       25396, 7390, 25486, 918,
                                                                       939, 8506, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27322, 0, 3,
                                                                       25486, 7435, 25576, 939,
                                                                       960, 8569, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27448, 0, 3,
                                                                       25576, 7480, 25666, 960,
                                                                       981, 8632, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27574, 0, 3,
                                                                       25666, 7525, 25756, 981,
                                                                       1002, 8695, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27700, 0, 3,
                                                                       25756, 7570, 25846, 1002,
                                                                       1023, 8758, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27826, 0, 3,
                                                                       25846, 7615, 25936, 1023,
                                                                       1044, 8821, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27952, 0, 3,
                                                                       25936, 7660, 26026, 1044,
                                                                       1065, 8884, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 28078, 0, 3,
                                                                       26026, 7705, 26116, 1065,
                                                                       1086, 8947, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 28204, 0, 3,
                                                                       26116, 7750, 26206, 1086,
                                                                       1107, 9010, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 28330, 0, 3,
                                                                       26296, 7930, 26386, 1149,
                                                                       1170, 9199, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 28456, 0, 3,
                                                                       26386, 7975, 26476, 1170,
                                                                       1191, 9262, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 28582, 0, 3,
                                                                       26476, 8020, 26566, 1191,
                                                                       1212, 9325, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 28708, 0, 3,
                                                                       26566, 8065, 26656, 1212,
                                                                       1233, 9388, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 28834, 0, 3,
                                                                       26656, 8110, 26746, 1233,
                                                                       1254, 9451, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 28960, 0, 3,
                                                                       26746, 8155, 26836, 1254,
                                                                       1275, 9514, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 29086, 0, 3,
                                                                       26836, 8200, 26926, 1275,
                                                                       1296, 9577, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 29212, 0, 3,
                                                                       26926, 8245, 27016, 1296,
                                                                       1317, 9640, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 29338, 0, 3,
                                                                       27016, 8290, 27106, 1317,
                                                                       1338, 9703, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 29464, 0, 3,
                                                                       27196, 8506, 27322, 1380,
                                                                       1408, 9934, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 29632, 0, 3,
                                                                       27322, 8569, 27448, 1408,
                                                                       1436, 10018, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 29800, 0, 3,
                                                                       27448, 8632, 27574, 1436,
                                                                       1464, 10102, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 29968, 0, 3,
                                                                       27574, 8695, 27700, 1464,
                                                                       1492, 10186, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 30136, 0, 3,
                                                                       27700, 8758, 27826, 1492,
                                                                       1520, 10270, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 30304, 0, 3,
                                                                       27826, 8821, 27952, 1520,
                                                                       1548, 10354, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 30472, 0, 3,
                                                                       27952, 8884, 28078, 1548,
                                                                       1576, 10438, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 30640, 0, 3,
                                                                       28078, 8947, 28204, 1576,
                                                                       1604, 10522, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 30808, 0, 3,
                                                                       28330, 9199, 28456, 1660,
                                                                       1688, 10774, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 30976, 0, 3,
                                                                       28456, 9262, 28582, 1688,
                                                                       1716, 10858, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 31144, 0, 3,
                                                                       28582, 9325, 28708, 1716,
                                                                       1744, 10942, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 31312, 0, 3,
                                                                       28708, 9388, 28834, 1744,
                                                                       1772, 11026, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 31480, 0, 3,
                                                                       28834, 9451, 28960, 1772,
                                                                       1800, 11110, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 31648, 0, 3,
                                                                       28960, 9514, 29086, 1800,
                                                                       1828, 11194, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 31816, 0, 3,
                                                                       29086, 9577, 29212, 1828,
                                                                       1856, 11278, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 31984, 0, 3,
                                                                       29212, 9640, 29338, 1856,
                                                                       1884, 11362, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 32152, 0, 3,
                                                                       29464, 9934, 29632, 1940,
                                                                       1976, 11662, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 32368, 0, 3,
                                                                       29632, 10018, 29800, 1976,
                                                                       2012, 11770, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 32584, 0, 3,
                                                                       29800, 10102, 29968, 2012,
                                                                       2048, 11878, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 32800, 0, 3,
                                                                       29968, 10186, 30136, 2048,
                                                                       2084, 11986, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 33016, 0, 3,
                                                                       30136, 10270, 30304, 2084,
                                                                       2120, 12094, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 33232, 0, 3,
                                                                       30304, 10354, 30472, 2120,
                                                                       2156, 12202, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 33448, 0, 3,
                                                                       30472, 10438, 30640, 2156,
                                                                       2192, 12310, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 33664, 0, 3,
                                                                       30808, 10774, 30976, 2264,
                                                                       2300, 12634, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 33880, 0, 3,
                                                                       30976, 10858, 31144, 2300,
                                                                       2336, 12742, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 34096, 0, 3,
                                                                       31144, 10942, 31312, 2336,
                                                                       2372, 12850, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 34312, 0, 3,
                                                                       31312, 11026, 31480, 2372,
                                                                       2408, 12958, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 34528, 0, 3,
                                                                       31480, 11110, 31648, 2408,
                                                                       2444, 13066, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 34744, 0, 3,
                                                                       31648, 11194, 31816, 2444,
                                                                       2480, 13174, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 34960, 0, 3,
                                                                       31816, 11278, 31984, 2480,
                                                                       2516, 13282, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 35176, 0, 3,
                                                                       32152, 11662, 32368, 2588,
                                                                       2633, 13660, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 35446, 0, 3,
                                                                       32368, 11770, 32584, 2633,
                                                                       2678, 13795, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 35716, 0, 3,
                                                                       32584, 11878, 32800, 2678,
                                                                       2723, 13930, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 35986, 0, 3,
                                                                       32800, 11986, 33016, 2723,
                                                                       2768, 14065, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 36256, 0, 3,
                                                                       33016, 12094, 33232, 2768,
                                                                       2813, 14200, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 36526, 0, 3,
                                                                       33232, 12202, 33448, 2813,
                                                                       2858, 14335, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 36796, 0, 3,
                                                                       33664, 12634, 33880, 2948,
                                                                       2993, 14740, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 37066, 0, 3,
                                                                       33880, 12742, 34096, 2993,
                                                                       3038, 14875, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 37336, 0, 3,
                                                                       34096, 12850, 34312, 3038,
                                                                       3083, 15010, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 37606, 0, 3,
                                                                       34312, 12958, 34528, 3083,
                                                                       3128, 15145, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 37876, 0, 3,
                                                                       34528, 13066, 34744, 3128,
                                                                       3173, 15280, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 38146, 0, 3,
                                                                       34744, 13174, 34960, 3173,
                                                                       3218, 15415, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 38416, 0, 3,
                                                                       35176, 13660, 35446, 3308,
                                                                       3363, 15880, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 38746, 0, 3,
                                                                       35446, 13795, 35716, 3363,
                                                                       3418, 16045, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 39076, 0, 3,
                                                                       35716, 13930, 35986, 3418,
                                                                       3473, 16210, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 39406, 0, 3,
                                                                       35986, 14065, 36256, 3473,
                                                                       3528, 16375, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 39736, 0, 3,
                                                                       36256, 14200, 36526, 3528,
                                                                       3583, 16540, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 40066, 0, 3,
                                                                       36796, 14740, 37066, 3693,
                                                                       3748, 17035, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 40396, 0, 3,
                                                                       37066, 14875, 37336, 3748,
                                                                       3803, 17200, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 40726, 0, 3,
                                                                       37336, 15010, 37606, 3803,
                                                                       3858, 17365, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 41056, 0, 3,
                                                                       37606, 15145, 37876, 3858,
                                                                       3913, 17530, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 41386, 0, 3,
                                                                       37876, 15280, 38146, 3913,
                                                                       3968, 17695, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 41716, 0, 3,
                                                                       38416, 15880, 38746, 4078,
                                                                       4144, 18256, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 42112, 0, 3,
                                                                       38746, 16045, 39076, 4144,
                                                                       4210, 18454, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 42508, 0, 3,
                                                                       39076, 16210, 39406, 4210,
                                                                       4276, 18652, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 42904, 0, 3,
                                                                       39406, 16375, 39736, 4276,
                                                                       4342, 18850, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 43300, 0, 3,
                                                                       40066, 17035, 40396, 4474,
                                                                       4540, 19444, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 43696, 0, 3,
                                                                       40396, 17200, 40726, 4540,
                                                                       4606, 19642, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 44092, 0, 3,
                                                                       40726, 17365, 41056, 4606,
                                                                       4672, 19840, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 44488, 0, 3,
                                                                       41056, 17530, 41386, 4672,
                                                                       4738, 20038, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 44884, 0, 3,
                                                                       41716, 18256, 42112, 4870,
                                                                       4948, 20704, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 45352, 0, 3,
                                                                       42112, 18454, 42508, 4948,
                                                                       5026, 20938, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 45820, 0, 3,
                                                                       42508, 18652, 42904, 5026,
                                                                       5104, 21172, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 46288, 0, 3,
                                                                       43300, 19444, 43696, 5260,
                                                                       5338, 21874, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 46756, 0, 3,
                                                                       43696, 19642, 44092, 5338,
                                                                       5416, 22108, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 47224, 0, 3,
                                                                       44092, 19840, 44488, 5416,
                                                                       5494, 22342, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47692, 3, 5650,
                                                                       5653, 22576, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47702, 3, 5653,
                                                                       5656, 22582, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47712, 3, 5656,
                                                                       5659, 22588, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47722, 3, 5659,
                                                                       5662, 22594, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47732, 3, 5662,
                                                                       5665, 22600, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47742, 3, 5665,
                                                                       5668, 22606, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47752, 3, 5668,
                                                                       5671, 22612, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47762, 3, 5671,
                                                                       5674, 22618, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47772, 3, 5674,
                                                                       5677, 22624, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47782, 3, 5677,
                                                                       5680, 22630, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47792, 3, 5680,
                                                                       5683, 22636, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47802, 3, 5683,
                                                                       5686, 22642, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47812, 3, 5686,
                                                                       5689, 22648, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47822, 3, 5689,
                                                                       5692, 22654, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47832, 3, 5698,
                                                                       5701, 22660, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47842, 3, 5701,
                                                                       5704, 22666, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47852, 3, 5704,
                                                                       5707, 22672, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47862, 3, 5707,
                                                                       5710, 22678, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47872, 3, 5710,
                                                                       5713, 22684, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47882, 3, 5713,
                                                                       5716, 22690, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47892, 3, 5716,
                                                                       5719, 22696, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47902, 3, 5719,
                                                                       5722, 22702, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47912, 3, 5722,
                                                                       5725, 22708, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47922, 3, 5725,
                                                                       5728, 22714, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47932, 3, 5728,
                                                                       5731, 22720, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47942, 3, 5731,
                                                                       5734, 22726, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47952, 3, 5734,
                                                                       5737, 22732, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47962, 3, 5737,
                                                                       5740, 22738, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 47972, 0, 3,
                                                                       47692, 22576, 47702, 5746,
                                                                       5755, 22744, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48002, 0, 3,
                                                                       47702, 22582, 47712, 5755,
                                                                       5764, 22762, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48032, 0, 3,
                                                                       47712, 22588, 47722, 5764,
                                                                       5773, 22780, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48062, 0, 3,
                                                                       47722, 22594, 47732, 5773,
                                                                       5782, 22798, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48092, 0, 3,
                                                                       47732, 22600, 47742, 5782,
                                                                       5791, 22816, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48122, 0, 3,
                                                                       47742, 22606, 47752, 5791,
                                                                       5800, 22834, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48152, 0, 3,
                                                                       47752, 22612, 47762, 5800,
                                                                       5809, 22852, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48182, 0, 3,
                                                                       47762, 22618, 47772, 5809,
                                                                       5818, 22870, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48212, 0, 3,
                                                                       47772, 22624, 47782, 5818,
                                                                       5827, 22888, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48242, 0, 3,
                                                                       47782, 22630, 47792, 5827,
                                                                       5836, 22906, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48272, 0, 3,
                                                                       47792, 22636, 47802, 5836,
                                                                       5845, 22924, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48302, 0, 3,
                                                                       47802, 22642, 47812, 5845,
                                                                       5854, 22942, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48332, 0, 3,
                                                                       47812, 22648, 47822, 5854,
                                                                       5863, 22960, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48362, 0, 3,
                                                                       47832, 22660, 47842, 5881,
                                                                       5890, 22978, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48392, 0, 3,
                                                                       47842, 22666, 47852, 5890,
                                                                       5899, 22996, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48422, 0, 3,
                                                                       47852, 22672, 47862, 5899,
                                                                       5908, 23014, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48452, 0, 3,
                                                                       47862, 22678, 47872, 5908,
                                                                       5917, 23032, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48482, 0, 3,
                                                                       47872, 22684, 47882, 5917,
                                                                       5926, 23050, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48512, 0, 3,
                                                                       47882, 22690, 47892, 5926,
                                                                       5935, 23068, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48542, 0, 3,
                                                                       47892, 22696, 47902, 5935,
                                                                       5944, 23086, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48572, 0, 3,
                                                                       47902, 22702, 47912, 5944,
                                                                       5953, 23104, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48602, 0, 3,
                                                                       47912, 22708, 47922, 5953,
                                                                       5962, 23122, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48632, 0, 3,
                                                                       47922, 22714, 47932, 5962,
                                                                       5971, 23140, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48662, 0, 3,
                                                                       47932, 22720, 47942, 5971,
                                                                       5980, 23158, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48692, 0, 3,
                                                                       47942, 22726, 47952, 5980,
                                                                       5989, 23176, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 48722, 0, 3,
                                                                       47952, 22732, 47962, 5989,
                                                                       5998, 23194, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 48752, 0, 3,
                                                                       47972, 22744, 48002, 6016,
                                                                       6034, 23212, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 48812, 0, 3,
                                                                       48002, 22762, 48032, 6034,
                                                                       6052, 23248, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 48872, 0, 3,
                                                                       48032, 22780, 48062, 6052,
                                                                       6070, 23284, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 48932, 0, 3,
                                                                       48062, 22798, 48092, 6070,
                                                                       6088, 23320, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 48992, 0, 3,
                                                                       48092, 22816, 48122, 6088,
                                                                       6106, 23356, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49052, 0, 3,
                                                                       48122, 22834, 48152, 6106,
                                                                       6124, 23392, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49112, 0, 3,
                                                                       48152, 22852, 48182, 6124,
                                                                       6142, 23428, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49172, 0, 3,
                                                                       48182, 22870, 48212, 6142,
                                                                       6160, 23464, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49232, 0, 3,
                                                                       48212, 22888, 48242, 6160,
                                                                       6178, 23500, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49292, 0, 3,
                                                                       48242, 22906, 48272, 6178,
                                                                       6196, 23536, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49352, 0, 3,
                                                                       48272, 22924, 48302, 6196,
                                                                       6214, 23572, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49412, 0, 3,
                                                                       48302, 22942, 48332, 6214,
                                                                       6232, 23608, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49472, 0, 3,
                                                                       48362, 22978, 48392, 6268,
                                                                       6286, 23644, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49532, 0, 3,
                                                                       48392, 22996, 48422, 6286,
                                                                       6304, 23680, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49592, 0, 3,
                                                                       48422, 23014, 48452, 6304,
                                                                       6322, 23716, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49652, 0, 3,
                                                                       48452, 23032, 48482, 6322,
                                                                       6340, 23752, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49712, 0, 3,
                                                                       48482, 23050, 48512, 6340,
                                                                       6358, 23788, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49772, 0, 3,
                                                                       48512, 23068, 48542, 6358,
                                                                       6376, 23824, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49832, 0, 3,
                                                                       48542, 23086, 48572, 6376,
                                                                       6394, 23860, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49892, 0, 3,
                                                                       48572, 23104, 48602, 6394,
                                                                       6412, 23896, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 49952, 0, 3,
                                                                       48602, 23122, 48632, 6412,
                                                                       6430, 23932, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 50012, 0, 3,
                                                                       48632, 23140, 48662, 6430,
                                                                       6448, 23968, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 50072, 0, 3,
                                                                       48662, 23158, 48692, 6448,
                                                                       6466, 24004, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 50132, 0, 3,
                                                                       48692, 23176, 48722, 6466,
                                                                       6484, 24040, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 50192, 0, 3,
                                                                       48752, 23212, 48812, 6520,
                                                                       6550, 24076, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 50292, 0, 3,
                                                                       48812, 23248, 48872, 6550,
                                                                       6580, 24136, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 50392, 0, 3,
                                                                       48872, 23284, 48932, 6580,
                                                                       6610, 24196, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 50492, 0, 3,
                                                                       48932, 23320, 48992, 6610,
                                                                       6640, 24256, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 50592, 0, 3,
                                                                       48992, 23356, 49052, 6640,
                                                                       6670, 24316, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 50692, 0, 3,
                                                                       49052, 23392, 49112, 6670,
                                                                       6700, 24376, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 50792, 0, 3,
                                                                       49112, 23428, 49172, 6700,
                                                                       6730, 24436, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 50892, 0, 3,
                                                                       49172, 23464, 49232, 6730,
                                                                       6760, 24496, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 50992, 0, 3,
                                                                       49232, 23500, 49292, 6760,
                                                                       6790, 24556, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 51092, 0, 3,
                                                                       49292, 23536, 49352, 6790,
                                                                       6820, 24616, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 51192, 0, 3,
                                                                       49352, 23572, 49412, 6820,
                                                                       6850, 24676, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 51292, 0, 3,
                                                                       49472, 23644, 49532, 6910,
                                                                       6940, 24736, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 51392, 0, 3,
                                                                       49532, 23680, 49592, 6940,
                                                                       6970, 24796, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 51492, 0, 3,
                                                                       49592, 23716, 49652, 6970,
                                                                       7000, 24856, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 51592, 0, 3,
                                                                       49652, 23752, 49712, 7000,
                                                                       7030, 24916, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 51692, 0, 3,
                                                                       49712, 23788, 49772, 7030,
                                                                       7060, 24976, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 51792, 0, 3,
                                                                       49772, 23824, 49832, 7060,
                                                                       7090, 25036, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 51892, 0, 3,
                                                                       49832, 23860, 49892, 7090,
                                                                       7120, 25096, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 51992, 0, 3,
                                                                       49892, 23896, 49952, 7120,
                                                                       7150, 25156, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 52092, 0, 3,
                                                                       49952, 23932, 50012, 7150,
                                                                       7180, 25216, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 52192, 0, 3,
                                                                       50012, 23968, 50072, 7180,
                                                                       7210, 25276, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 52292, 0, 3,
                                                                       50072, 24004, 50132, 7210,
                                                                       7240, 25336, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 52392, 0, 3,
                                                                       50192, 24076, 50292, 7300,
                                                                       7345, 25396, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 52542, 0, 3,
                                                                       50292, 24136, 50392, 7345,
                                                                       7390, 25486, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 52692, 0, 3,
                                                                       50392, 24196, 50492, 7390,
                                                                       7435, 25576, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 52842, 0, 3,
                                                                       50492, 24256, 50592, 7435,
                                                                       7480, 25666, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 52992, 0, 3,
                                                                       50592, 24316, 50692, 7480,
                                                                       7525, 25756, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 53142, 0, 3,
                                                                       50692, 24376, 50792, 7525,
                                                                       7570, 25846, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 53292, 0, 3,
                                                                       50792, 24436, 50892, 7570,
                                                                       7615, 25936, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 53442, 0, 3,
                                                                       50892, 24496, 50992, 7615,
                                                                       7660, 26026, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 53592, 0, 3,
                                                                       50992, 24556, 51092, 7660,
                                                                       7705, 26116, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 53742, 0, 3,
                                                                       51092, 24616, 51192, 7705,
                                                                       7750, 26206, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 53892, 0, 3,
                                                                       51292, 24736, 51392, 7840,
                                                                       7885, 26296, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 54042, 0, 3,
                                                                       51392, 24796, 51492, 7885,
                                                                       7930, 26386, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 54192, 0, 3,
                                                                       51492, 24856, 51592, 7930,
                                                                       7975, 26476, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 54342, 0, 3,
                                                                       51592, 24916, 51692, 7975,
                                                                       8020, 26566, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 54492, 0, 3,
                                                                       51692, 24976, 51792, 8020,
                                                                       8065, 26656, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 54642, 0, 3,
                                                                       51792, 25036, 51892, 8065,
                                                                       8110, 26746, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 54792, 0, 3,
                                                                       51892, 25096, 51992, 8110,
                                                                       8155, 26836, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 54942, 0, 3,
                                                                       51992, 25156, 52092, 8155,
                                                                       8200, 26926, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 55092, 0, 3,
                                                                       52092, 25216, 52192, 8200,
                                                                       8245, 27016, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 55242, 0, 3,
                                                                       52192, 25276, 52292, 8245,
                                                                       8290, 27106, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 55392, 0, 3,
                                                                       52392, 25396, 52542, 8380,
                                                                       8443, 27196, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 55602, 0, 3,
                                                                       52542, 25486, 52692, 8443,
                                                                       8506, 27322, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 55812, 0, 3,
                                                                       52692, 25576, 52842, 8506,
                                                                       8569, 27448, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 56022, 0, 3,
                                                                       52842, 25666, 52992, 8569,
                                                                       8632, 27574, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 56232, 0, 3,
                                                                       52992, 25756, 53142, 8632,
                                                                       8695, 27700, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 56442, 0, 3,
                                                                       53142, 25846, 53292, 8695,
                                                                       8758, 27826, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 56652, 0, 3,
                                                                       53292, 25936, 53442, 8758,
                                                                       8821, 27952, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 56862, 0, 3,
                                                                       53442, 26026, 53592, 8821,
                                                                       8884, 28078, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 57072, 0, 3,
                                                                       53592, 26116, 53742, 8884,
                                                                       8947, 28204, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 57282, 0, 3,
                                                                       53892, 26296, 54042, 9073,
                                                                       9136, 28330, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 57492, 0, 3,
                                                                       54042, 26386, 54192, 9136,
                                                                       9199, 28456, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 57702, 0, 3,
                                                                       54192, 26476, 54342, 9199,
                                                                       9262, 28582, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 57912, 0, 3,
                                                                       54342, 26566, 54492, 9262,
                                                                       9325, 28708, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 58122, 0, 3,
                                                                       54492, 26656, 54642, 9325,
                                                                       9388, 28834, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 58332, 0, 3,
                                                                       54642, 26746, 54792, 9388,
                                                                       9451, 28960, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 58542, 0, 3,
                                                                       54792, 26836, 54942, 9451,
                                                                       9514, 29086, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 58752, 0, 3,
                                                                       54942, 26926, 55092, 9514,
                                                                       9577, 29212, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 58962, 0, 3,
                                                                       55092, 27016, 55242, 9577,
                                                                       9640, 29338, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 59172, 0, 3,
                                                                       55392, 27196, 55602, 9766,
                                                                       9850, 29464, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 59452, 0, 3,
                                                                       55602, 27322, 55812, 9850,
                                                                       9934, 29632, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 59732, 0, 3,
                                                                       55812, 27448, 56022, 9934,
                                                                       10018, 29800, ncols,
                                                                       gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 60012, 0, 3,
                                                                       56022, 27574, 56232,
                                                                       10018, 10102, 29968,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 60292, 0, 3,
                                                                       56232, 27700, 56442,
                                                                       10102, 10186, 30136,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 60572, 0, 3,
                                                                       56442, 27826, 56652,
                                                                       10186, 10270, 30304,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 60852, 0, 3,
                                                                       56652, 27952, 56862,
                                                                       10270, 10354, 30472,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 61132, 0, 3,
                                                                       56862, 28078, 57072,
                                                                       10354, 10438, 30640,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 61412, 0, 3,
                                                                       57282, 28330, 57492,
                                                                       10606, 10690, 30808,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 61692, 0, 3,
                                                                       57492, 28456, 57702,
                                                                       10690, 10774, 30976,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 61972, 0, 3,
                                                                       57702, 28582, 57912,
                                                                       10774, 10858, 31144,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 62252, 0, 3,
                                                                       57912, 28708, 58122,
                                                                       10858, 10942, 31312,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 62532, 0, 3,
                                                                       58122, 28834, 58332,
                                                                       10942, 11026, 31480,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 62812, 0, 3,
                                                                       58332, 28960, 58542,
                                                                       11026, 11110, 31648,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 63092, 0, 3,
                                                                       58542, 29086, 58752,
                                                                       11110, 11194, 31816,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 63372, 0, 3,
                                                                       58752, 29212, 58962,
                                                                       11194, 11278, 31984,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 63652, 0, 3,
                                                                       59172, 29464, 59452,
                                                                       11446, 11554, 32152,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 64012, 0, 3,
                                                                       59452, 29632, 59732,
                                                                       11554, 11662, 32368,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 64372, 0, 3,
                                                                       59732, 29800, 60012,
                                                                       11662, 11770, 32584,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 64732, 0, 3,
                                                                       60012, 29968, 60292,
                                                                       11770, 11878, 32800,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 65092, 0, 3,
                                                                       60292, 30136, 60572,
                                                                       11878, 11986, 33016,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 65452, 0, 3,
                                                                       60572, 30304, 60852,
                                                                       11986, 12094, 33232,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 65812, 0, 3,
                                                                       60852, 30472, 61132,
                                                                       12094, 12202, 33448,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 66172, 0, 3,
                                                                       61412, 30808, 61692,
                                                                       12418, 12526, 33664,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 66532, 0, 3,
                                                                       61692, 30976, 61972,
                                                                       12526, 12634, 33880,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 66892, 0, 3,
                                                                       61972, 31144, 62252,
                                                                       12634, 12742, 34096,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 67252, 0, 3,
                                                                       62252, 31312, 62532,
                                                                       12742, 12850, 34312,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 67612, 0, 3,
                                                                       62532, 31480, 62812,
                                                                       12850, 12958, 34528,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 67972, 0, 3,
                                                                       62812, 31648, 63092,
                                                                       12958, 13066, 34744,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 68332, 0, 3,
                                                                       63092, 31816, 63372,
                                                                       13066, 13174, 34960,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 68692, 0, 3,
                                                                       63652, 32152, 64012,
                                                                       13390, 13525, 35176,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 69142, 0, 3,
                                                                       64012, 32368, 64372,
                                                                       13525, 13660, 35446,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 69592, 0, 3,
                                                                       64372, 32584, 64732,
                                                                       13660, 13795, 35716,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 70042, 0, 3,
                                                                       64732, 32800, 65092,
                                                                       13795, 13930, 35986,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 70492, 0, 3,
                                                                       65092, 33016, 65452,
                                                                       13930, 14065, 36256,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 70942, 0, 3,
                                                                       65452, 33232, 65812,
                                                                       14065, 14200, 36526,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 71392, 0, 3,
                                                                       66172, 33664, 66532,
                                                                       14470, 14605, 36796,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 71842, 0, 3,
                                                                       66532, 33880, 66892,
                                                                       14605, 14740, 37066,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 72292, 0, 3,
                                                                       66892, 34096, 67252,
                                                                       14740, 14875, 37336,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 72742, 0, 3,
                                                                       67252, 34312, 67612,
                                                                       14875, 15010, 37606,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 73192, 0, 3,
                                                                       67612, 34528, 67972,
                                                                       15010, 15145, 37876,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 73642, 0, 3,
                                                                       67972, 34744, 68332,
                                                                       15145, 15280, 38146,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 74092, 0, 3,
                                                                       68692, 35176, 69142,
                                                                       15550, 15715, 38416,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 74642, 0, 3,
                                                                       69142, 35446, 69592,
                                                                       15715, 15880, 38746,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 75192, 0, 3,
                                                                       69592, 35716, 70042,
                                                                       15880, 16045, 39076,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 75742, 0, 3,
                                                                       70042, 35986, 70492,
                                                                       16045, 16210, 39406,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 76292, 0, 3,
                                                                       70492, 36256, 70942,
                                                                       16210, 16375, 39736,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 76842, 0, 3,
                                                                       71392, 36796, 71842,
                                                                       16705, 16870, 40066,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 77392, 0, 3,
                                                                       71842, 37066, 72292,
                                                                       16870, 17035, 40396,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 77942, 0, 3,
                                                                       72292, 37336, 72742,
                                                                       17035, 17200, 40726,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 78492, 0, 3,
                                                                       72742, 37606, 73192,
                                                                       17200, 17365, 41056,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 79042, 0, 3,
                                                                       73192, 37876, 73642,
                                                                       17365, 17530, 41386,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 79592, 0, 3,
                                                                       74092, 38416, 74642,
                                                                       17860, 18058, 41716,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 80252, 0, 3,
                                                                       74642, 38746, 75192,
                                                                       18058, 18256, 42112,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 80912, 0, 3,
                                                                       75192, 39076, 75742,
                                                                       18256, 18454, 42508,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 81572, 0, 3,
                                                                       75742, 39406, 76292,
                                                                       18454, 18652, 42904,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 82232, 0, 3,
                                                                       76842, 40066, 77392,
                                                                       19048, 19246, 43300,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 82892, 0, 3,
                                                                       77392, 40396, 77942,
                                                                       19246, 19444, 43696,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 83552, 0, 3,
                                                                       77942, 40726, 78492,
                                                                       19444, 19642, 44092,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 84212, 0, 3,
                                                                       78492, 41056, 79042,
                                                                       19642, 19840, 44488,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 84872, 0, 3,
                                                                       79592, 41716, 80252,
                                                                       20236, 20470, 44884,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 85652, 0, 3,
                                                                       80252, 42112, 80912,
                                                                       20470, 20704, 45352,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 86432, 0, 3,
                                                                       80912, 42508, 81572,
                                                                       20704, 20938, 45820,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 87212, 0, 3,
                                                                       82232, 43300, 82892,
                                                                       21406, 21640, 46288,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 87992, 0, 3,
                                                                       82892, 43696, 83552,
                                                                       21640, 21874, 46756,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 88772, 0, 3,
                                                                       83552, 44092, 84212,
                                                                       21874, 22108, 47224,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89552, 3, 22576,
                                                                       22582, 47712, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89567, 3, 22582,
                                                                       22588, 47722, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89582, 3, 22588,
                                                                       22594, 47732, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89597, 3, 22594,
                                                                       22600, 47742, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89612, 3, 22600,
                                                                       22606, 47752, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89627, 3, 22606,
                                                                       22612, 47762, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89642, 3, 22612,
                                                                       22618, 47772, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89657, 3, 22618,
                                                                       22624, 47782, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89672, 3, 22624,
                                                                       22630, 47792, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89687, 3, 22630,
                                                                       22636, 47802, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89702, 3, 22636,
                                                                       22642, 47812, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89717, 3, 22642,
                                                                       22648, 47822, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89732, 3, 22660,
                                                                       22666, 47852, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89747, 3, 22666,
                                                                       22672, 47862, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89762, 3, 22672,
                                                                       22678, 47872, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89777, 3, 22678,
                                                                       22684, 47882, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89792, 3, 22684,
                                                                       22690, 47892, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89807, 3, 22690,
                                                                       22696, 47902, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89822, 3, 22696,
                                                                       22702, 47912, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89837, 3, 22702,
                                                                       22708, 47922, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89852, 3, 22708,
                                                                       22714, 47932, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89867, 3, 22714,
                                                                       22720, 47942, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89882, 3, 22720,
                                                                       22726, 47952, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89897, 3, 22726,
                                                                       22732, 47962, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 89912, 0, 3,
                                                                       89552, 47712, 89567,
                                                                       22744, 22762, 48032,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 89957, 0, 3,
                                                                       89567, 47722, 89582,
                                                                       22762, 22780, 48062,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90002, 0, 3,
                                                                       89582, 47732, 89597,
                                                                       22780, 22798, 48092,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90047, 0, 3,
                                                                       89597, 47742, 89612,
                                                                       22798, 22816, 48122,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90092, 0, 3,
                                                                       89612, 47752, 89627,
                                                                       22816, 22834, 48152,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90137, 0, 3,
                                                                       89627, 47762, 89642,
                                                                       22834, 22852, 48182,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90182, 0, 3,
                                                                       89642, 47772, 89657,
                                                                       22852, 22870, 48212,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90227, 0, 3,
                                                                       89657, 47782, 89672,
                                                                       22870, 22888, 48242,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90272, 0, 3,
                                                                       89672, 47792, 89687,
                                                                       22888, 22906, 48272,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90317, 0, 3,
                                                                       89687, 47802, 89702,
                                                                       22906, 22924, 48302,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90362, 0, 3,
                                                                       89702, 47812, 89717,
                                                                       22924, 22942, 48332,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90407, 0, 3,
                                                                       89732, 47852, 89747,
                                                                       22978, 22996, 48422,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90452, 0, 3,
                                                                       89747, 47862, 89762,
                                                                       22996, 23014, 48452,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90497, 0, 3,
                                                                       89762, 47872, 89777,
                                                                       23014, 23032, 48482,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90542, 0, 3,
                                                                       89777, 47882, 89792,
                                                                       23032, 23050, 48512,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90587, 0, 3,
                                                                       89792, 47892, 89807,
                                                                       23050, 23068, 48542,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90632, 0, 3,
                                                                       89807, 47902, 89822,
                                                                       23068, 23086, 48572,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90677, 0, 3,
                                                                       89822, 47912, 89837,
                                                                       23086, 23104, 48602,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90722, 0, 3,
                                                                       89837, 47922, 89852,
                                                                       23104, 23122, 48632,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90767, 0, 3,
                                                                       89852, 47932, 89867,
                                                                       23122, 23140, 48662,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90812, 0, 3,
                                                                       89867, 47942, 89882,
                                                                       23140, 23158, 48692,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 90857, 0, 3,
                                                                       89882, 47952, 89897,
                                                                       23158, 23176, 48722,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 90902, 0, 3,
                                                                       89912, 48032, 89957,
                                                                       23212, 23248, 48872,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 90992, 0, 3,
                                                                       89957, 48062, 90002,
                                                                       23248, 23284, 48932,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 91082, 0, 3,
                                                                       90002, 48092, 90047,
                                                                       23284, 23320, 48992,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 91172, 0, 3,
                                                                       90047, 48122, 90092,
                                                                       23320, 23356, 49052,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 91262, 0, 3,
                                                                       90092, 48152, 90137,
                                                                       23356, 23392, 49112,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 91352, 0, 3,
                                                                       90137, 48182, 90182,
                                                                       23392, 23428, 49172,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 91442, 0, 3,
                                                                       90182, 48212, 90227,
                                                                       23428, 23464, 49232,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 91532, 0, 3,
                                                                       90227, 48242, 90272,
                                                                       23464, 23500, 49292,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 91622, 0, 3,
                                                                       90272, 48272, 90317,
                                                                       23500, 23536, 49352,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 91712, 0, 3,
                                                                       90317, 48302, 90362,
                                                                       23536, 23572, 49412,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 91802, 0, 3,
                                                                       90407, 48422, 90452,
                                                                       23644, 23680, 49592,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 91892, 0, 3,
                                                                       90452, 48452, 90497,
                                                                       23680, 23716, 49652,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 91982, 0, 3,
                                                                       90497, 48482, 90542,
                                                                       23716, 23752, 49712,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 92072, 0, 3,
                                                                       90542, 48512, 90587,
                                                                       23752, 23788, 49772,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 92162, 0, 3,
                                                                       90587, 48542, 90632,
                                                                       23788, 23824, 49832,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 92252, 0, 3,
                                                                       90632, 48572, 90677,
                                                                       23824, 23860, 49892,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 92342, 0, 3,
                                                                       90677, 48602, 90722,
                                                                       23860, 23896, 49952,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 92432, 0, 3,
                                                                       90722, 48632, 90767,
                                                                       23896, 23932, 50012,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 92522, 0, 3,
                                                                       90767, 48662, 90812,
                                                                       23932, 23968, 50072,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 92612, 0, 3,
                                                                       90812, 48692, 90857,
                                                                       23968, 24004, 50132,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 92702, 0, 3,
                                                                       90902, 48872, 90992,
                                                                       24076, 24136, 50392,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 92852, 0, 3,
                                                                       90992, 48932, 91082,
                                                                       24136, 24196, 50492,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 93002, 0, 3,
                                                                       91082, 48992, 91172,
                                                                       24196, 24256, 50592,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 93152, 0, 3,
                                                                       91172, 49052, 91262,
                                                                       24256, 24316, 50692,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 93302, 0, 3,
                                                                       91262, 49112, 91352,
                                                                       24316, 24376, 50792,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 93452, 0, 3,
                                                                       91352, 49172, 91442,
                                                                       24376, 24436, 50892,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 93602, 0, 3,
                                                                       91442, 49232, 91532,
                                                                       24436, 24496, 50992,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 93752, 0, 3,
                                                                       91532, 49292, 91622,
                                                                       24496, 24556, 51092,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 93902, 0, 3,
                                                                       91622, 49352, 91712,
                                                                       24556, 24616, 51192,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 94052, 0, 3,
                                                                       91802, 49592, 91892,
                                                                       24736, 24796, 51492,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 94202, 0, 3,
                                                                       91892, 49652, 91982,
                                                                       24796, 24856, 51592,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 94352, 0, 3,
                                                                       91982, 49712, 92072,
                                                                       24856, 24916, 51692,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 94502, 0, 3,
                                                                       92072, 49772, 92162,
                                                                       24916, 24976, 51792,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 94652, 0, 3,
                                                                       92162, 49832, 92252,
                                                                       24976, 25036, 51892,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 94802, 0, 3,
                                                                       92252, 49892, 92342,
                                                                       25036, 25096, 51992,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 94952, 0, 3,
                                                                       92342, 49952, 92432,
                                                                       25096, 25156, 52092,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 95102, 0, 3,
                                                                       92432, 50012, 92522,
                                                                       25156, 25216, 52192,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 95252, 0, 3,
                                                                       92522, 50072, 92612,
                                                                       25216, 25276, 52292,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 95402, 0, 3,
                                                                       92702, 50392, 92852,
                                                                       25396, 25486, 52692,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 95627, 0, 3,
                                                                       92852, 50492, 93002,
                                                                       25486, 25576, 52842,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 95852, 0, 3,
                                                                       93002, 50592, 93152,
                                                                       25576, 25666, 52992,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 96077, 0, 3,
                                                                       93152, 50692, 93302,
                                                                       25666, 25756, 53142,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 96302, 0, 3,
                                                                       93302, 50792, 93452,
                                                                       25756, 25846, 53292,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 96527, 0, 3,
                                                                       93452, 50892, 93602,
                                                                       25846, 25936, 53442,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 96752, 0, 3,
                                                                       93602, 50992, 93752,
                                                                       25936, 26026, 53592,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 96977, 0, 3,
                                                                       93752, 51092, 93902,
                                                                       26026, 26116, 53742,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 97202, 0, 3,
                                                                       94052, 51492, 94202,
                                                                       26296, 26386, 54192,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 97427, 0, 3,
                                                                       94202, 51592, 94352,
                                                                       26386, 26476, 54342,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 97652, 0, 3,
                                                                       94352, 51692, 94502,
                                                                       26476, 26566, 54492,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 97877, 0, 3,
                                                                       94502, 51792, 94652,
                                                                       26566, 26656, 54642,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 98102, 0, 3,
                                                                       94652, 51892, 94802,
                                                                       26656, 26746, 54792,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 98327, 0, 3,
                                                                       94802, 51992, 94952,
                                                                       26746, 26836, 54942,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 98552, 0, 3,
                                                                       94952, 52092, 95102,
                                                                       26836, 26926, 55092,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 98777, 0, 3,
                                                                       95102, 52192, 95252,
                                                                       26926, 27016, 55242,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 99002, 0, 3,
                                                                       95402, 52692, 95627,
                                                                       27196, 27322, 55812,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 99317, 0, 3,
                                                                       95627, 52842, 95852,
                                                                       27322, 27448, 56022,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 99632, 0, 3,
                                                                       95852, 52992, 96077,
                                                                       27448, 27574, 56232,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 99947, 0, 3,
                                                                       96077, 53142, 96302,
                                                                       27574, 27700, 56442,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 100262, 0, 3,
                                                                       96302, 53292, 96527,
                                                                       27700, 27826, 56652,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 100577, 0, 3,
                                                                       96527, 53442, 96752,
                                                                       27826, 27952, 56862,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 100892, 0, 3,
                                                                       96752, 53592, 96977,
                                                                       27952, 28078, 57072,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 101207, 0, 3,
                                                                       97202, 54192, 97427,
                                                                       28330, 28456, 57702,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 101522, 0, 3,
                                                                       97427, 54342, 97652,
                                                                       28456, 28582, 57912,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 101837, 0, 3,
                                                                       97652, 54492, 97877,
                                                                       28582, 28708, 58122,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 102152, 0, 3,
                                                                       97877, 54642, 98102,
                                                                       28708, 28834, 58332,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 102467, 0, 3,
                                                                       98102, 54792, 98327,
                                                                       28834, 28960, 58542,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 102782, 0, 3,
                                                                       98327, 54942, 98552,
                                                                       28960, 29086, 58752,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 103097, 0, 3,
                                                                       98552, 55092, 98777,
                                                                       29086, 29212, 58962,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 103412, 0, 3,
                                                                       99002, 55812, 99317,
                                                                       29464, 29632, 59732,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 103832, 0, 3,
                                                                       99317, 56022, 99632,
                                                                       29632, 29800, 60012,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 104252, 0, 3,
                                                                       99632, 56232, 99947,
                                                                       29800, 29968, 60292,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 104672, 0, 3,
                                                                       99947, 56442, 100262,
                                                                       29968, 30136, 60572,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 105092, 0, 3,
                                                                       100262, 56652, 100577,
                                                                       30136, 30304, 60852,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 105512, 0, 3,
                                                                       100577, 56862, 100892,
                                                                       30304, 30472, 61132,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 105932, 0, 3,
                                                                       101207, 57702, 101522,
                                                                       30808, 30976, 61972,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 106352, 0, 3,
                                                                       101522, 57912, 101837,
                                                                       30976, 31144, 62252,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 106772, 0, 3,
                                                                       101837, 58122, 102152,
                                                                       31144, 31312, 62532,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 107192, 0, 3,
                                                                       102152, 58332, 102467,
                                                                       31312, 31480, 62812,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 107612, 0, 3,
                                                                       102467, 58542, 102782,
                                                                       31480, 31648, 63092,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 108032, 0, 3,
                                                                       102782, 58752, 103097,
                                                                       31648, 31816, 63372,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 108452, 0, 3,
                                                                       103412, 59732, 103832,
                                                                       32152, 32368, 64372,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 108992, 0, 3,
                                                                       103832, 60012, 104252,
                                                                       32368, 32584, 64732,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 109532, 0, 3,
                                                                       104252, 60292, 104672,
                                                                       32584, 32800, 65092,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 110072, 0, 3,
                                                                       104672, 60572, 105092,
                                                                       32800, 33016, 65452,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 110612, 0, 3,
                                                                       105092, 60852, 105512,
                                                                       33016, 33232, 65812,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 111152, 0, 3,
                                                                       105932, 61972, 106352,
                                                                       33664, 33880, 66892,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 111692, 0, 3,
                                                                       106352, 62252, 106772,
                                                                       33880, 34096, 67252,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 112232, 0, 3,
                                                                       106772, 62532, 107192,
                                                                       34096, 34312, 67612,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 112772, 0, 3,
                                                                       107192, 62812, 107612,
                                                                       34312, 34528, 67972,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 113312, 0, 3,
                                                                       107612, 63092, 108032,
                                                                       34528, 34744, 68332,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 113852, 0, 3,
                                                                       108452, 64372, 108992,
                                                                       35176, 35446, 69592,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 114527, 0, 3,
                                                                       108992, 64732, 109532,
                                                                       35446, 35716, 70042,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 115202, 0, 3,
                                                                       109532, 65092, 110072,
                                                                       35716, 35986, 70492,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 115877, 0, 3,
                                                                       110072, 65452, 110612,
                                                                       35986, 36256, 70942,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 116552, 0, 3,
                                                                       111152, 66892, 111692,
                                                                       36796, 37066, 72292,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 117227, 0, 3,
                                                                       111692, 67252, 112232,
                                                                       37066, 37336, 72742,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 117902, 0, 3,
                                                                       112232, 67612, 112772,
                                                                       37336, 37606, 73192,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 118577, 0, 3,
                                                                       112772, 67972, 113312,
                                                                       37606, 37876, 73642,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 119252, 0, 3,
                                                                       113852, 69592, 114527,
                                                                       38416, 38746, 75192,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 120077, 0, 3,
                                                                       114527, 70042, 115202,
                                                                       38746, 39076, 75742,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 120902, 0, 3,
                                                                       115202, 70492, 115877,
                                                                       39076, 39406, 76292,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 121727, 0, 3,
                                                                       116552, 72292, 117227,
                                                                       40066, 40396, 77942,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 122552, 0, 3,
                                                                       117227, 72742, 117902,
                                                                       40396, 40726, 78492,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 123377, 0, 3,
                                                                       117902, 73192, 118577,
                                                                       40726, 41056, 79042,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 124202, 0, 3,
                                                                       119252, 75192, 120077,
                                                                       41716, 42112, 80912,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 125192, 0, 3,
                                                                       120077, 75742, 120902,
                                                                       42112, 42508, 81572,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 126182, 0, 3,
                                                                       121727, 77942, 122552,
                                                                       43300, 43696, 83552,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 127172, 0, 3,
                                                                       122552, 78492, 123377,
                                                                       43696, 44092, 84212,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 128162, 0, 3,
                                                                       124202, 80912, 125192,
                                                                       44884, 45352, 86432,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 129332, 0, 3,
                                                                       126182, 83552, 127172,
                                                                       46288, 46756, 88772,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130502, 3, 47692,
                                                                       47702, 89552, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130523, 3, 47702,
                                                                       47712, 89567, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130544, 3, 47712,
                                                                       47722, 89582, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130565, 3, 47722,
                                                                       47732, 89597, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130586, 3, 47732,
                                                                       47742, 89612, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130607, 3, 47742,
                                                                       47752, 89627, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130628, 3, 47752,
                                                                       47762, 89642, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130649, 3, 47762,
                                                                       47772, 89657, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130670, 3, 47772,
                                                                       47782, 89672, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130691, 3, 47782,
                                                                       47792, 89687, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130712, 3, 47792,
                                                                       47802, 89702, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130733, 3, 47802,
                                                                       47812, 89717, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130754, 3, 47832,
                                                                       47842, 89732, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130775, 3, 47842,
                                                                       47852, 89747, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130796, 3, 47852,
                                                                       47862, 89762, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130817, 3, 47862,
                                                                       47872, 89777, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130838, 3, 47872,
                                                                       47882, 89792, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130859, 3, 47882,
                                                                       47892, 89807, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130880, 3, 47892,
                                                                       47902, 89822, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130901, 3, 47902,
                                                                       47912, 89837, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130922, 3, 47912,
                                                                       47922, 89852, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130943, 3, 47922,
                                                                       47932, 89867, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130964, 3, 47932,
                                                                       47942, 89882, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130985, 3, 47942,
                                                                       47952, 89897, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 131006, 0, 3,
                                                                       130502, 89552, 130523,
                                                                       47972, 48002, 89912,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 131069, 0, 3,
                                                                       130523, 89567, 130544,
                                                                       48002, 48032, 89957,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 131132, 0, 3,
                                                                       130544, 89582, 130565,
                                                                       48032, 48062, 90002,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 131195, 0, 3,
                                                                       130565, 89597, 130586,
                                                                       48062, 48092, 90047,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 131258, 0, 3,
                                                                       130586, 89612, 130607,
                                                                       48092, 48122, 90092,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 131321, 0, 3,
                                                                       130607, 89627, 130628,
                                                                       48122, 48152, 90137,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 131384, 0, 3,
                                                                       130628, 89642, 130649,
                                                                       48152, 48182, 90182,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 131447, 0, 3,
                                                                       130649, 89657, 130670,
                                                                       48182, 48212, 90227,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 131510, 0, 3,
                                                                       130670, 89672, 130691,
                                                                       48212, 48242, 90272,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 131573, 0, 3,
                                                                       130691, 89687, 130712,
                                                                       48242, 48272, 90317,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 131636, 0, 3,
                                                                       130712, 89702, 130733,
                                                                       48272, 48302, 90362,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 131699, 0, 3,
                                                                       130754, 89732, 130775,
                                                                       48362, 48392, 90407,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 131762, 0, 3,
                                                                       130775, 89747, 130796,
                                                                       48392, 48422, 90452,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 131825, 0, 3,
                                                                       130796, 89762, 130817,
                                                                       48422, 48452, 90497,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 131888, 0, 3,
                                                                       130817, 89777, 130838,
                                                                       48452, 48482, 90542,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 131951, 0, 3,
                                                                       130838, 89792, 130859,
                                                                       48482, 48512, 90587,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 132014, 0, 3,
                                                                       130859, 89807, 130880,
                                                                       48512, 48542, 90632,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 132077, 0, 3,
                                                                       130880, 89822, 130901,
                                                                       48542, 48572, 90677,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 132140, 0, 3,
                                                                       130901, 89837, 130922,
                                                                       48572, 48602, 90722,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 132203, 0, 3,
                                                                       130922, 89852, 130943,
                                                                       48602, 48632, 90767,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 132266, 0, 3,
                                                                       130943, 89867, 130964,
                                                                       48632, 48662, 90812,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 132329, 0, 3,
                                                                       130964, 89882, 130985,
                                                                       48662, 48692, 90857,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 132392, 0, 3,
                                                                       131006, 89912, 131069,
                                                                       48752, 48812, 90902,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 132518, 0, 3,
                                                                       131069, 89957, 131132,
                                                                       48812, 48872, 90992,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 132644, 0, 3,
                                                                       131132, 90002, 131195,
                                                                       48872, 48932, 91082,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 132770, 0, 3,
                                                                       131195, 90047, 131258,
                                                                       48932, 48992, 91172,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 132896, 0, 3,
                                                                       131258, 90092, 131321,
                                                                       48992, 49052, 91262,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 133022, 0, 3,
                                                                       131321, 90137, 131384,
                                                                       49052, 49112, 91352,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 133148, 0, 3,
                                                                       131384, 90182, 131447,
                                                                       49112, 49172, 91442,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 133274, 0, 3,
                                                                       131447, 90227, 131510,
                                                                       49172, 49232, 91532,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 133400, 0, 3,
                                                                       131510, 90272, 131573,
                                                                       49232, 49292, 91622,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 133526, 0, 3,
                                                                       131573, 90317, 131636,
                                                                       49292, 49352, 91712,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 133652, 0, 3,
                                                                       131699, 90407, 131762,
                                                                       49472, 49532, 91802,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 133778, 0, 3,
                                                                       131762, 90452, 131825,
                                                                       49532, 49592, 91892,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 133904, 0, 3,
                                                                       131825, 90497, 131888,
                                                                       49592, 49652, 91982,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 134030, 0, 3,
                                                                       131888, 90542, 131951,
                                                                       49652, 49712, 92072,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 134156, 0, 3,
                                                                       131951, 90587, 132014,
                                                                       49712, 49772, 92162,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 134282, 0, 3,
                                                                       132014, 90632, 132077,
                                                                       49772, 49832, 92252,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 134408, 0, 3,
                                                                       132077, 90677, 132140,
                                                                       49832, 49892, 92342,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 134534, 0, 3,
                                                                       132140, 90722, 132203,
                                                                       49892, 49952, 92432,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 134660, 0, 3,
                                                                       132203, 90767, 132266,
                                                                       49952, 50012, 92522,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 134786, 0, 3,
                                                                       132266, 90812, 132329,
                                                                       50012, 50072, 92612,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 134912, 0, 3,
                                                                       132392, 90902, 132518,
                                                                       50192, 50292, 92702,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 135122, 0, 3,
                                                                       132518, 90992, 132644,
                                                                       50292, 50392, 92852,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 135332, 0, 3,
                                                                       132644, 91082, 132770,
                                                                       50392, 50492, 93002,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 135542, 0, 3,
                                                                       132770, 91172, 132896,
                                                                       50492, 50592, 93152,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 135752, 0, 3,
                                                                       132896, 91262, 133022,
                                                                       50592, 50692, 93302,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 135962, 0, 3,
                                                                       133022, 91352, 133148,
                                                                       50692, 50792, 93452,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 136172, 0, 3,
                                                                       133148, 91442, 133274,
                                                                       50792, 50892, 93602,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 136382, 0, 3,
                                                                       133274, 91532, 133400,
                                                                       50892, 50992, 93752,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 136592, 0, 3,
                                                                       133400, 91622, 133526,
                                                                       50992, 51092, 93902,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 136802, 0, 3,
                                                                       133652, 91802, 133778,
                                                                       51292, 51392, 94052,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 137012, 0, 3,
                                                                       133778, 91892, 133904,
                                                                       51392, 51492, 94202,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 137222, 0, 3,
                                                                       133904, 91982, 134030,
                                                                       51492, 51592, 94352,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 137432, 0, 3,
                                                                       134030, 92072, 134156,
                                                                       51592, 51692, 94502,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 137642, 0, 3,
                                                                       134156, 92162, 134282,
                                                                       51692, 51792, 94652,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 137852, 0, 3,
                                                                       134282, 92252, 134408,
                                                                       51792, 51892, 94802,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 138062, 0, 3,
                                                                       134408, 92342, 134534,
                                                                       51892, 51992, 94952,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 138272, 0, 3,
                                                                       134534, 92432, 134660,
                                                                       51992, 52092, 95102,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 138482, 0, 3,
                                                                       134660, 92522, 134786,
                                                                       52092, 52192, 95252,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 138692, 0, 3,
                                                                       134912, 92702, 135122,
                                                                       52392, 52542, 95402,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 139007, 0, 3,
                                                                       135122, 92852, 135332,
                                                                       52542, 52692, 95627,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 139322, 0, 3,
                                                                       135332, 93002, 135542,
                                                                       52692, 52842, 95852,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 139637, 0, 3,
                                                                       135542, 93152, 135752,
                                                                       52842, 52992, 96077,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 139952, 0, 3,
                                                                       135752, 93302, 135962,
                                                                       52992, 53142, 96302,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 140267, 0, 3,
                                                                       135962, 93452, 136172,
                                                                       53142, 53292, 96527,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 140582, 0, 3,
                                                                       136172, 93602, 136382,
                                                                       53292, 53442, 96752,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 140897, 0, 3,
                                                                       136382, 93752, 136592,
                                                                       53442, 53592, 96977,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 141212, 0, 3,
                                                                       136802, 94052, 137012,
                                                                       53892, 54042, 97202,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 141527, 0, 3,
                                                                       137012, 94202, 137222,
                                                                       54042, 54192, 97427,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 141842, 0, 3,
                                                                       137222, 94352, 137432,
                                                                       54192, 54342, 97652,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 142157, 0, 3,
                                                                       137432, 94502, 137642,
                                                                       54342, 54492, 97877,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 142472, 0, 3,
                                                                       137642, 94652, 137852,
                                                                       54492, 54642, 98102,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 142787, 0, 3,
                                                                       137852, 94802, 138062,
                                                                       54642, 54792, 98327,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 143102, 0, 3,
                                                                       138062, 94952, 138272,
                                                                       54792, 54942, 98552,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 143417, 0, 3,
                                                                       138272, 95102, 138482,
                                                                       54942, 55092, 98777,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 143732, 0, 3,
                                                                       138692, 95402, 139007,
                                                                       55392, 55602, 99002,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 144173, 0, 3,
                                                                       139007, 95627, 139322,
                                                                       55602, 55812, 99317,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 144614, 0, 3,
                                                                       139322, 95852, 139637,
                                                                       55812, 56022, 99632,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 145055, 0, 3,
                                                                       139637, 96077, 139952,
                                                                       56022, 56232, 99947,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 145496, 0, 3,
                                                                       139952, 96302, 140267,
                                                                       56232, 56442, 100262,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 145937, 0, 3,
                                                                       140267, 96527, 140582,
                                                                       56442, 56652, 100577,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 146378, 0, 3,
                                                                       140582, 96752, 140897,
                                                                       56652, 56862, 100892,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 146819, 0, 3,
                                                                       141212, 97202, 141527,
                                                                       57282, 57492, 101207,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 147260, 0, 3,
                                                                       141527, 97427, 141842,
                                                                       57492, 57702, 101522,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 147701, 0, 3,
                                                                       141842, 97652, 142157,
                                                                       57702, 57912, 101837,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 148142, 0, 3,
                                                                       142157, 97877, 142472,
                                                                       57912, 58122, 102152,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 148583, 0, 3,
                                                                       142472, 98102, 142787,
                                                                       58122, 58332, 102467,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 149024, 0, 3,
                                                                       142787, 98327, 143102,
                                                                       58332, 58542, 102782,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 149465, 0, 3,
                                                                       143102, 98552, 143417,
                                                                       58542, 58752, 103097,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 149906, 0, 3,
                                                                       143732, 99002, 144173,
                                                                       59172, 59452, 103412,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 150494, 0, 3,
                                                                       144173, 99317, 144614,
                                                                       59452, 59732, 103832,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 151082, 0, 3,
                                                                       144614, 99632, 145055,
                                                                       59732, 60012, 104252,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 151670, 0, 3,
                                                                       145055, 99947, 145496,
                                                                       60012, 60292, 104672,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 152258, 0, 3,
                                                                       145496, 100262, 145937,
                                                                       60292, 60572, 105092,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 152846, 0, 3,
                                                                       145937, 100577, 146378,
                                                                       60572, 60852, 105512,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 153434, 0, 3,
                                                                       146819, 101207, 147260,
                                                                       61412, 61692, 105932,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 154022, 0, 3,
                                                                       147260, 101522, 147701,
                                                                       61692, 61972, 106352,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 154610, 0, 3,
                                                                       147701, 101837, 148142,
                                                                       61972, 62252, 106772,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 155198, 0, 3,
                                                                       148142, 102152, 148583,
                                                                       62252, 62532, 107192,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 155786, 0, 3,
                                                                       148583, 102467, 149024,
                                                                       62532, 62812, 107612,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 156374, 0, 3,
                                                                       149024, 102782, 149465,
                                                                       62812, 63092, 108032,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 156962, 0, 3,
                                                                       149906, 103412, 150494,
                                                                       63652, 64012, 108452,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 157718, 0, 3,
                                                                       150494, 103832, 151082,
                                                                       64012, 64372, 108992,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 158474, 0, 3,
                                                                       151082, 104252, 151670,
                                                                       64372, 64732, 109532,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 159230, 0, 3,
                                                                       151670, 104672, 152258,
                                                                       64732, 65092, 110072,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 159986, 0, 3,
                                                                       152258, 105092, 152846,
                                                                       65092, 65452, 110612,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 160742, 0, 3,
                                                                       153434, 105932, 154022,
                                                                       66172, 66532, 111152,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 161498, 0, 3,
                                                                       154022, 106352, 154610,
                                                                       66532, 66892, 111692,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 162254, 0, 3,
                                                                       154610, 106772, 155198,
                                                                       66892, 67252, 112232,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 163010, 0, 3,
                                                                       155198, 107192, 155786,
                                                                       67252, 67612, 112772,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 163766, 0, 3,
                                                                       155786, 107612, 156374,
                                                                       67612, 67972, 113312,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 164522, 0, 3,
                                                                       156962, 108452, 157718,
                                                                       68692, 69142, 113852,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 165467, 0, 3,
                                                                       157718, 108992, 158474,
                                                                       69142, 69592, 114527,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 166412, 0, 3,
                                                                       158474, 109532, 159230,
                                                                       69592, 70042, 115202,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 167357, 0, 3,
                                                                       159230, 110072, 159986,
                                                                       70042, 70492, 115877,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 168302, 0, 3,
                                                                       160742, 111152, 161498,
                                                                       71392, 71842, 116552,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 169247, 0, 3,
                                                                       161498, 111692, 162254,
                                                                       71842, 72292, 117227,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 170192, 0, 3,
                                                                       162254, 112232, 163010,
                                                                       72292, 72742, 117902,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 171137, 0, 3,
                                                                       163010, 112772, 163766,
                                                                       72742, 73192, 118577,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 172082, 0, 3,
                                                                       164522, 113852, 165467,
                                                                       74092, 74642, 119252,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 173237, 0, 3,
                                                                       165467, 114527, 166412,
                                                                       74642, 75192, 120077,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 174392, 0, 3,
                                                                       166412, 115202, 167357,
                                                                       75192, 75742, 120902,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 175547, 0, 3,
                                                                       168302, 116552, 169247,
                                                                       76842, 77392, 121727,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 176702, 0, 3,
                                                                       169247, 117227, 170192,
                                                                       77392, 77942, 122552,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 177857, 0, 3,
                                                                       170192, 117902, 171137,
                                                                       77942, 78492, 123377,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 179012, 0, 3,
                                                                       172082, 119252, 173237,
                                                                       79592, 80252, 124202,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 180398, 0, 3,
                                                                       173237, 120077, 174392,
                                                                       80252, 80912, 125192,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 181784, 0, 3,
                                                                       175547, 121727, 176702,
                                                                       82232, 82892, 126182,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 183170, 0, 3,
                                                                       176702, 122552, 177857,
                                                                       82892, 83552, 127172,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 184556, 0, 3,
                                                                       179012, 124202, 180398,
                                                                       84872, 85652, 128162,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 186194, 0, 3,
                                                                       181784, 126182, 183170,
                                                                       87212, 87992, 129332,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 187832, 149906, 588, ncols);

                    simdfunc::contract_primitives(buffer, 188728, 153434, 588, ncols);

                    simdfunc::contract_primitives(buffer, 189624, 156962, 756, ncols);

                    simdfunc::contract_primitives(buffer, 190776, 160742, 756, ncols);

                    simdfunc::contract_primitives(buffer, 191928, 164522, 945, ncols);

                    simdfunc::contract_primitives(buffer, 193368, 168302, 945, ncols);

                    simdfunc::contract_primitives(buffer, 194808, 172082, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 196568, 175547, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 198328, 179012, 1386, ncols);

                    simdfunc::contract_primitives(buffer, 200440, 181784, 1386, ncols);

                    simdfunc::contract_primitives(buffer, 202552, 184556, 1638, ncols);

                    simdfunc::contract_primitives(buffer, 205048, 186194, 1638, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 188420, 187832, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 189316, 188728, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 190380, 189624, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 191532, 190776, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 192873, 191928, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 194313, 193368, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 195963, 194808, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 197723, 196568, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 199714, 198328, 66, 1, nmax);

        simdtrf::transform_h_inner(buffer, 201826, 200440, 66, 1, nmax);

        simdtrf::transform_h_inner(buffer, 204190, 202552, 78, 1, nmax);

        simdtrf::transform_h_inner(buffer, 206686, 205048, 78, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 207544, 188420, 190380, 11,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 208468, 189316, 191532, 11,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 209392, 190380, 192873, 11,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 210580, 191532, 194313, 11,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 211768, 192873, 195963, 11,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 213253, 194313, 197723, 11,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 214738, 195963, 199714, 11,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 216553, 197723, 201826, 11,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 218368, 199714, 204190, 11,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 220546, 201826, 206686, 11,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 222724, 207544, 209392, 11,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 224572, 208468, 210580, 11,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 226420, 209392, 211768, 11,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 228796, 210580, 213253, 11,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 231172, 211768, 214738, 11,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 234142, 213253, 216553, 11,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 237112, 214738, 218368, 11,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 240742, 216553, 220546, 11,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 244372, 222724, 226420, 11,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 247452, 224572, 228796, 11,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 250532, 226420, 231172, 11,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 254492, 228796, 234142, 11,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 258452, 231172, 237112, 11,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 263402, 234142, 240742, 11,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 268352, 244372, 250532, 11,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 272972, 247452, 254492, 11,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 277592, 250532, 258452, 11,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 283532, 254492, 263402, 11,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 289472, 268352, 277592, 11,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 295940, 272972, 283532, 11,
                                             nmax);

        simdtrf::transform_h_inner(buffer, 302408, 295940, 28, 11, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 302408, 121, nmax);

        simdtrf::transform_h_inner(buffer, 302408, 289472, 28, 11, nmax);

        simdtrf::transform_i_outer(values + 1573 * nvalues + n * npairs, nvalues, buffer, 302408,
                                   121, nmax);
    }

    for (size_t m = 0; m < 3146; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
