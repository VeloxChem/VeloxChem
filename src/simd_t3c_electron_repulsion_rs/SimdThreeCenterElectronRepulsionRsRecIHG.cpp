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


#include "SimdThreeCenterElectronRepulsionRsRecIHG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
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
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_ihg_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_ihg_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 206746, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2574 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 206746, 111574, 14082, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 15,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 23, 3, 15,
                                                             ncols, fj, i * nprim_b + j, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5650, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5653, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5656, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5659, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5662, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5665, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5668, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5671, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5674, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5677, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5680, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5683, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5686, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5689, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5692, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5695, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5698, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5701, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5704, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5707, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5710, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5713, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5716, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5719, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5722, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5725, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5728, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5731, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5734, 3, 9, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5743, 3, 10, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5752, 3, 11, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5761, 3, 12, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5770, 3, 13, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5779, 3, 14, 61,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5788, 3, 15, 64,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5797, 3, 16, 67,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5806, 3, 17, 70,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5815, 3, 18, 73,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5824, 3, 19, 76,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5833, 3, 20, 79,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5842, 3, 21, 82,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5851, 3, 26, 91,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5860, 3, 27, 94,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5869, 3, 28, 97,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5878, 3, 29, 100,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5887, 3, 30, 103,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5896, 3, 31, 106,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5905, 3, 32, 109,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5914, 3, 33, 112,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5923, 3, 34, 115,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5932, 3, 35, 118,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5941, 3, 36, 121,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5950, 3, 37, 124,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5959, 3, 38, 127,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5968, 3, 46, 142,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5986, 3, 49, 148,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6004, 3, 52, 154,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6022, 3, 55, 160,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6040, 3, 58, 166,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6058, 3, 61, 172,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6076, 3, 64, 178,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6094, 3, 67, 184,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6112, 3, 70, 190,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6130, 3, 73, 196,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6148, 3, 76, 202,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6166, 3, 79, 208,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6184, 3, 91, 226,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6202, 3, 94, 232,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6220, 3, 97, 238,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6238, 3, 100, 244,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6256, 3, 103, 250,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6274, 3, 106, 256,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6292, 3, 109, 262,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6310, 3, 112, 268,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6328, 3, 115, 274,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6346, 3, 118, 280,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6364, 3, 121, 286,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6382, 3, 124, 292,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6400, 3, 142, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6430, 3, 148, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6460, 3, 154, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6490, 3, 160, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6520, 3, 166, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6550, 3, 172, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6580, 3, 178, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6610, 3, 184, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6640, 3, 190, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6670, 3, 196, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6700, 3, 202, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6730, 3, 226, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6760, 3, 232, 458,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6790, 3, 238, 468,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6820, 3, 244, 478,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6850, 3, 250, 488,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6880, 3, 256, 498,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6910, 3, 262, 508,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6940, 3, 268, 518,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6970, 3, 274, 528,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 7000, 3, 280, 538,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 7030, 3, 286, 548,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7060, 3, 318, 588,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7105, 3, 328, 603,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7150, 3, 338, 618,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7195, 3, 348, 633,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7240, 3, 358, 648,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7285, 3, 368, 663,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7330, 3, 378, 678,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7375, 3, 388, 693,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7420, 3, 398, 708,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7465, 3, 408, 723,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7510, 3, 448, 768,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7555, 3, 458, 783,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7600, 3, 468, 798,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7645, 3, 478, 813,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7690, 3, 488, 828,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7735, 3, 498, 843,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7780, 3, 508, 858,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7825, 3, 518, 873,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7870, 3, 528, 888,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7915, 3, 538, 903,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7960, 3, 588, 960,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8023, 3, 603, 981,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8086, 3, 618,
                                                                       1002, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8149, 3, 633,
                                                                       1023, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8212, 3, 648,
                                                                       1044, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8275, 3, 663,
                                                                       1065, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8338, 3, 678,
                                                                       1086, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8401, 3, 693,
                                                                       1107, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8464, 3, 708,
                                                                       1128, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8527, 3, 768,
                                                                       1191, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8590, 3, 783,
                                                                       1212, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8653, 3, 798,
                                                                       1233, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8716, 3, 813,
                                                                       1254, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8779, 3, 828,
                                                                       1275, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8842, 3, 843,
                                                                       1296, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8905, 3, 858,
                                                                       1317, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8968, 3, 873,
                                                                       1338, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9031, 3, 888,
                                                                       1359, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9094, 3, 960,
                                                                       1436, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9178, 3, 981,
                                                                       1464, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9262, 3, 1002,
                                                                       1492, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9346, 3, 1023,
                                                                       1520, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9430, 3, 1044,
                                                                       1548, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9514, 3, 1065,
                                                                       1576, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9598, 3, 1086,
                                                                       1604, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9682, 3, 1107,
                                                                       1632, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9766, 3, 1191,
                                                                       1716, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9850, 3, 1212,
                                                                       1744, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9934, 3, 1233,
                                                                       1772, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10018, 3, 1254,
                                                                       1800, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10102, 3, 1275,
                                                                       1828, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10186, 3, 1296,
                                                                       1856, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10270, 3, 1317,
                                                                       1884, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10354, 3, 1338,
                                                                       1912, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10438, 3, 1436,
                                                                       2012, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10546, 3, 1464,
                                                                       2048, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10654, 3, 1492,
                                                                       2084, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10762, 3, 1520,
                                                                       2120, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10870, 3, 1548,
                                                                       2156, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10978, 3, 1576,
                                                                       2192, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11086, 3, 1604,
                                                                       2228, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11194, 3, 1716,
                                                                       2336, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11302, 3, 1744,
                                                                       2372, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11410, 3, 1772,
                                                                       2408, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11518, 3, 1800,
                                                                       2444, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11626, 3, 1828,
                                                                       2480, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11734, 3, 1856,
                                                                       2516, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11842, 3, 1884,
                                                                       2552, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11950, 3, 2012,
                                                                       2678, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12085, 3, 2048,
                                                                       2723, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12220, 3, 2084,
                                                                       2768, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12355, 3, 2120,
                                                                       2813, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12490, 3, 2156,
                                                                       2858, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12625, 3, 2192,
                                                                       2903, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12760, 3, 2336,
                                                                       3038, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12895, 3, 2372,
                                                                       3083, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13030, 3, 2408,
                                                                       3128, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13165, 3, 2444,
                                                                       3173, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13300, 3, 2480,
                                                                       3218, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13435, 3, 2516,
                                                                       3263, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13570, 3, 2678,
                                                                       3418, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13735, 3, 2723,
                                                                       3473, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13900, 3, 2768,
                                                                       3528, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14065, 3, 2813,
                                                                       3583, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14230, 3, 2858,
                                                                       3638, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14395, 3, 3038,
                                                                       3803, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14560, 3, 3083,
                                                                       3858, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14725, 3, 3128,
                                                                       3913, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14890, 3, 3173,
                                                                       3968, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15055, 3, 3218,
                                                                       4023, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 15220, 3, 3418,
                                                                       4210, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 15418, 3, 3473,
                                                                       4276, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 15616, 3, 3528,
                                                                       4342, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 15814, 3, 3583,
                                                                       4408, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 16012, 3, 3803,
                                                                       4606, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 16210, 3, 3858,
                                                                       4672, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 16408, 3, 3913,
                                                                       4738, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 16606, 3, 3968,
                                                                       4804, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 16804, 3, 4210,
                                                                       5026, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 17038, 3, 4276,
                                                                       5104, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 17272, 3, 4342,
                                                                       5182, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 17506, 3, 4606,
                                                                       5416, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 17740, 3, 4672,
                                                                       5494, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 17974, 3, 4738,
                                                                       5572, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18208, 3, 7, 8,
                                                                       5650, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18214, 3, 8, 9,
                                                                       5653, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18220, 3, 9, 10,
                                                                       5656, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18226, 3, 10, 11,
                                                                       5659, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18232, 3, 11, 12,
                                                                       5662, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18238, 3, 12, 13,
                                                                       5665, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18244, 3, 13, 14,
                                                                       5668, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18250, 3, 14, 15,
                                                                       5671, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18256, 3, 15, 16,
                                                                       5674, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18262, 3, 16, 17,
                                                                       5677, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18268, 3, 17, 18,
                                                                       5680, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18274, 3, 18, 19,
                                                                       5683, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18280, 3, 19, 20,
                                                                       5686, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18286, 3, 20, 21,
                                                                       5689, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18292, 3, 24, 25,
                                                                       5692, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18298, 3, 25, 26,
                                                                       5695, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18304, 3, 26, 27,
                                                                       5698, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18310, 3, 27, 28,
                                                                       5701, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18316, 3, 28, 29,
                                                                       5704, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18322, 3, 29, 30,
                                                                       5707, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18328, 3, 30, 31,
                                                                       5710, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18334, 3, 31, 32,
                                                                       5713, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18340, 3, 32, 33,
                                                                       5716, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18346, 3, 33, 34,
                                                                       5719, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18352, 3, 34, 35,
                                                                       5722, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18358, 3, 35, 36,
                                                                       5725, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18364, 3, 36, 37,
                                                                       5728, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18370, 3, 37, 38,
                                                                       5731, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18376, 0, 3,
                                                                       18208, 5650, 18214, 5734,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18394, 0, 3,
                                                                       18214, 5653, 18220, 5743,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18412, 0, 3,
                                                                       18220, 5656, 18226, 5752,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18430, 0, 3,
                                                                       18226, 5659, 18232, 5761,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18448, 0, 3,
                                                                       18232, 5662, 18238, 5770,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18466, 0, 3,
                                                                       18238, 5665, 18244, 5779,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18484, 0, 3,
                                                                       18244, 5668, 18250, 5788,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18502, 0, 3,
                                                                       18250, 5671, 18256, 5797,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18520, 0, 3,
                                                                       18256, 5674, 18262, 5806,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18538, 0, 3,
                                                                       18262, 5677, 18268, 5815,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18556, 0, 3,
                                                                       18268, 5680, 18274, 5824,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18574, 0, 3,
                                                                       18274, 5683, 18280, 5833,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18592, 0, 3,
                                                                       18280, 5686, 18286, 5842,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18610, 0, 3,
                                                                       18292, 5692, 18298, 5851,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18628, 0, 3,
                                                                       18298, 5695, 18304, 5860,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18646, 0, 3,
                                                                       18304, 5698, 18310, 5869,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18664, 0, 3,
                                                                       18310, 5701, 18316, 5878,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18682, 0, 3,
                                                                       18316, 5704, 18322, 5887,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18700, 0, 3,
                                                                       18322, 5707, 18328, 5896,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18718, 0, 3,
                                                                       18328, 5710, 18334, 5905,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18736, 0, 3,
                                                                       18334, 5713, 18340, 5914,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18754, 0, 3,
                                                                       18340, 5716, 18346, 5923,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18772, 0, 3,
                                                                       18346, 5719, 18352, 5932,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18790, 0, 3,
                                                                       18352, 5722, 18358, 5941,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18808, 0, 3,
                                                                       18358, 5725, 18364, 5950,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18826, 0, 3,
                                                                       18364, 5728, 18370, 5959,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18844, 0, 3,
                                                                       18376, 5734, 18394, 130,
                                                                       136, 5968, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18880, 0, 3,
                                                                       18394, 5743, 18412, 136,
                                                                       142, 5986, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18916, 0, 3,
                                                                       18412, 5752, 18430, 142,
                                                                       148, 6004, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18952, 0, 3,
                                                                       18430, 5761, 18448, 148,
                                                                       154, 6022, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18988, 0, 3,
                                                                       18448, 5770, 18466, 154,
                                                                       160, 6040, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19024, 0, 3,
                                                                       18466, 5779, 18484, 160,
                                                                       166, 6058, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19060, 0, 3,
                                                                       18484, 5788, 18502, 166,
                                                                       172, 6076, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19096, 0, 3,
                                                                       18502, 5797, 18520, 172,
                                                                       178, 6094, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19132, 0, 3,
                                                                       18520, 5806, 18538, 178,
                                                                       184, 6112, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19168, 0, 3,
                                                                       18538, 5815, 18556, 184,
                                                                       190, 6130, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19204, 0, 3,
                                                                       18556, 5824, 18574, 190,
                                                                       196, 6148, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19240, 0, 3,
                                                                       18574, 5833, 18592, 196,
                                                                       202, 6166, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19276, 0, 3,
                                                                       18610, 5851, 18628, 214,
                                                                       220, 6184, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19312, 0, 3,
                                                                       18628, 5860, 18646, 220,
                                                                       226, 6202, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19348, 0, 3,
                                                                       18646, 5869, 18664, 226,
                                                                       232, 6220, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19384, 0, 3,
                                                                       18664, 5878, 18682, 232,
                                                                       238, 6238, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19420, 0, 3,
                                                                       18682, 5887, 18700, 238,
                                                                       244, 6256, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19456, 0, 3,
                                                                       18700, 5896, 18718, 244,
                                                                       250, 6274, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19492, 0, 3,
                                                                       18718, 5905, 18736, 250,
                                                                       256, 6292, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19528, 0, 3,
                                                                       18736, 5914, 18754, 256,
                                                                       262, 6310, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19564, 0, 3,
                                                                       18754, 5923, 18772, 262,
                                                                       268, 6328, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19600, 0, 3,
                                                                       18772, 5932, 18790, 268,
                                                                       274, 6346, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19636, 0, 3,
                                                                       18790, 5941, 18808, 274,
                                                                       280, 6364, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19672, 0, 3,
                                                                       18808, 5950, 18826, 280,
                                                                       286, 6382, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19708, 0, 3,
                                                                       18844, 5968, 18880, 298,
                                                                       308, 6400, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19768, 0, 3,
                                                                       18880, 5986, 18916, 308,
                                                                       318, 6430, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19828, 0, 3,
                                                                       18916, 6004, 18952, 318,
                                                                       328, 6460, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19888, 0, 3,
                                                                       18952, 6022, 18988, 328,
                                                                       338, 6490, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19948, 0, 3,
                                                                       18988, 6040, 19024, 338,
                                                                       348, 6520, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20008, 0, 3,
                                                                       19024, 6058, 19060, 348,
                                                                       358, 6550, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20068, 0, 3,
                                                                       19060, 6076, 19096, 358,
                                                                       368, 6580, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20128, 0, 3,
                                                                       19096, 6094, 19132, 368,
                                                                       378, 6610, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20188, 0, 3,
                                                                       19132, 6112, 19168, 378,
                                                                       388, 6640, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20248, 0, 3,
                                                                       19168, 6130, 19204, 388,
                                                                       398, 6670, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20308, 0, 3,
                                                                       19204, 6148, 19240, 398,
                                                                       408, 6700, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20368, 0, 3,
                                                                       19276, 6184, 19312, 428,
                                                                       438, 6730, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20428, 0, 3,
                                                                       19312, 6202, 19348, 438,
                                                                       448, 6760, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20488, 0, 3,
                                                                       19348, 6220, 19384, 448,
                                                                       458, 6790, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20548, 0, 3,
                                                                       19384, 6238, 19420, 458,
                                                                       468, 6820, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20608, 0, 3,
                                                                       19420, 6256, 19456, 468,
                                                                       478, 6850, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20668, 0, 3,
                                                                       19456, 6274, 19492, 478,
                                                                       488, 6880, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20728, 0, 3,
                                                                       19492, 6292, 19528, 488,
                                                                       498, 6910, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20788, 0, 3,
                                                                       19528, 6310, 19564, 498,
                                                                       508, 6940, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20848, 0, 3,
                                                                       19564, 6328, 19600, 508,
                                                                       518, 6970, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20908, 0, 3,
                                                                       19600, 6346, 19636, 518,
                                                                       528, 7000, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20968, 0, 3,
                                                                       19636, 6364, 19672, 528,
                                                                       538, 7030, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21028, 0, 3,
                                                                       19708, 6400, 19768, 558,
                                                                       573, 7060, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21118, 0, 3,
                                                                       19768, 6430, 19828, 573,
                                                                       588, 7105, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21208, 0, 3,
                                                                       19828, 6460, 19888, 588,
                                                                       603, 7150, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21298, 0, 3,
                                                                       19888, 6490, 19948, 603,
                                                                       618, 7195, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21388, 0, 3,
                                                                       19948, 6520, 20008, 618,
                                                                       633, 7240, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21478, 0, 3,
                                                                       20008, 6550, 20068, 633,
                                                                       648, 7285, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21568, 0, 3,
                                                                       20068, 6580, 20128, 648,
                                                                       663, 7330, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21658, 0, 3,
                                                                       20128, 6610, 20188, 663,
                                                                       678, 7375, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21748, 0, 3,
                                                                       20188, 6640, 20248, 678,
                                                                       693, 7420, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21838, 0, 3,
                                                                       20248, 6670, 20308, 693,
                                                                       708, 7465, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21928, 0, 3,
                                                                       20368, 6730, 20428, 738,
                                                                       753, 7510, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22018, 0, 3,
                                                                       20428, 6760, 20488, 753,
                                                                       768, 7555, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22108, 0, 3,
                                                                       20488, 6790, 20548, 768,
                                                                       783, 7600, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22198, 0, 3,
                                                                       20548, 6820, 20608, 783,
                                                                       798, 7645, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22288, 0, 3,
                                                                       20608, 6850, 20668, 798,
                                                                       813, 7690, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22378, 0, 3,
                                                                       20668, 6880, 20728, 813,
                                                                       828, 7735, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22468, 0, 3,
                                                                       20728, 6910, 20788, 828,
                                                                       843, 7780, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22558, 0, 3,
                                                                       20788, 6940, 20848, 843,
                                                                       858, 7825, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22648, 0, 3,
                                                                       20848, 6970, 20908, 858,
                                                                       873, 7870, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22738, 0, 3,
                                                                       20908, 7000, 20968, 873,
                                                                       888, 7915, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22828, 0, 3,
                                                                       21028, 7060, 21118, 918,
                                                                       939, 7960, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22954, 0, 3,
                                                                       21118, 7105, 21208, 939,
                                                                       960, 8023, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23080, 0, 3,
                                                                       21208, 7150, 21298, 960,
                                                                       981, 8086, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23206, 0, 3,
                                                                       21298, 7195, 21388, 981,
                                                                       1002, 8149, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23332, 0, 3,
                                                                       21388, 7240, 21478, 1002,
                                                                       1023, 8212, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23458, 0, 3,
                                                                       21478, 7285, 21568, 1023,
                                                                       1044, 8275, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23584, 0, 3,
                                                                       21568, 7330, 21658, 1044,
                                                                       1065, 8338, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23710, 0, 3,
                                                                       21658, 7375, 21748, 1065,
                                                                       1086, 8401, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23836, 0, 3,
                                                                       21748, 7420, 21838, 1086,
                                                                       1107, 8464, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23962, 0, 3,
                                                                       21928, 7510, 22018, 1149,
                                                                       1170, 8527, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24088, 0, 3,
                                                                       22018, 7555, 22108, 1170,
                                                                       1191, 8590, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24214, 0, 3,
                                                                       22108, 7600, 22198, 1191,
                                                                       1212, 8653, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24340, 0, 3,
                                                                       22198, 7645, 22288, 1212,
                                                                       1233, 8716, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24466, 0, 3,
                                                                       22288, 7690, 22378, 1233,
                                                                       1254, 8779, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24592, 0, 3,
                                                                       22378, 7735, 22468, 1254,
                                                                       1275, 8842, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24718, 0, 3,
                                                                       22468, 7780, 22558, 1275,
                                                                       1296, 8905, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24844, 0, 3,
                                                                       22558, 7825, 22648, 1296,
                                                                       1317, 8968, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24970, 0, 3,
                                                                       22648, 7870, 22738, 1317,
                                                                       1338, 9031, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25096, 0, 3,
                                                                       22828, 7960, 22954, 1380,
                                                                       1408, 9094, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25264, 0, 3,
                                                                       22954, 8023, 23080, 1408,
                                                                       1436, 9178, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25432, 0, 3,
                                                                       23080, 8086, 23206, 1436,
                                                                       1464, 9262, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25600, 0, 3,
                                                                       23206, 8149, 23332, 1464,
                                                                       1492, 9346, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25768, 0, 3,
                                                                       23332, 8212, 23458, 1492,
                                                                       1520, 9430, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25936, 0, 3,
                                                                       23458, 8275, 23584, 1520,
                                                                       1548, 9514, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26104, 0, 3,
                                                                       23584, 8338, 23710, 1548,
                                                                       1576, 9598, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26272, 0, 3,
                                                                       23710, 8401, 23836, 1576,
                                                                       1604, 9682, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26440, 0, 3,
                                                                       23962, 8527, 24088, 1660,
                                                                       1688, 9766, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26608, 0, 3,
                                                                       24088, 8590, 24214, 1688,
                                                                       1716, 9850, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26776, 0, 3,
                                                                       24214, 8653, 24340, 1716,
                                                                       1744, 9934, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26944, 0, 3,
                                                                       24340, 8716, 24466, 1744,
                                                                       1772, 10018, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 27112, 0, 3,
                                                                       24466, 8779, 24592, 1772,
                                                                       1800, 10102, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 27280, 0, 3,
                                                                       24592, 8842, 24718, 1800,
                                                                       1828, 10186, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 27448, 0, 3,
                                                                       24718, 8905, 24844, 1828,
                                                                       1856, 10270, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 27616, 0, 3,
                                                                       24844, 8968, 24970, 1856,
                                                                       1884, 10354, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27784, 0, 3,
                                                                       25096, 9094, 25264, 1940,
                                                                       1976, 10438, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28000, 0, 3,
                                                                       25264, 9178, 25432, 1976,
                                                                       2012, 10546, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28216, 0, 3,
                                                                       25432, 9262, 25600, 2012,
                                                                       2048, 10654, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28432, 0, 3,
                                                                       25600, 9346, 25768, 2048,
                                                                       2084, 10762, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28648, 0, 3,
                                                                       25768, 9430, 25936, 2084,
                                                                       2120, 10870, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28864, 0, 3,
                                                                       25936, 9514, 26104, 2120,
                                                                       2156, 10978, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 29080, 0, 3,
                                                                       26104, 9598, 26272, 2156,
                                                                       2192, 11086, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 29296, 0, 3,
                                                                       26440, 9766, 26608, 2264,
                                                                       2300, 11194, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 29512, 0, 3,
                                                                       26608, 9850, 26776, 2300,
                                                                       2336, 11302, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 29728, 0, 3,
                                                                       26776, 9934, 26944, 2336,
                                                                       2372, 11410, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 29944, 0, 3,
                                                                       26944, 10018, 27112, 2372,
                                                                       2408, 11518, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 30160, 0, 3,
                                                                       27112, 10102, 27280, 2408,
                                                                       2444, 11626, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 30376, 0, 3,
                                                                       27280, 10186, 27448, 2444,
                                                                       2480, 11734, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 30592, 0, 3,
                                                                       27448, 10270, 27616, 2480,
                                                                       2516, 11842, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 30808, 0, 3,
                                                                       27784, 10438, 28000, 2588,
                                                                       2633, 11950, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 31078, 0, 3,
                                                                       28000, 10546, 28216, 2633,
                                                                       2678, 12085, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 31348, 0, 3,
                                                                       28216, 10654, 28432, 2678,
                                                                       2723, 12220, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 31618, 0, 3,
                                                                       28432, 10762, 28648, 2723,
                                                                       2768, 12355, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 31888, 0, 3,
                                                                       28648, 10870, 28864, 2768,
                                                                       2813, 12490, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 32158, 0, 3,
                                                                       28864, 10978, 29080, 2813,
                                                                       2858, 12625, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 32428, 0, 3,
                                                                       29296, 11194, 29512, 2948,
                                                                       2993, 12760, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 32698, 0, 3,
                                                                       29512, 11302, 29728, 2993,
                                                                       3038, 12895, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 32968, 0, 3,
                                                                       29728, 11410, 29944, 3038,
                                                                       3083, 13030, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 33238, 0, 3,
                                                                       29944, 11518, 30160, 3083,
                                                                       3128, 13165, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 33508, 0, 3,
                                                                       30160, 11626, 30376, 3128,
                                                                       3173, 13300, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 33778, 0, 3,
                                                                       30376, 11734, 30592, 3173,
                                                                       3218, 13435, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 34048, 0, 3,
                                                                       30808, 11950, 31078, 3308,
                                                                       3363, 13570, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 34378, 0, 3,
                                                                       31078, 12085, 31348, 3363,
                                                                       3418, 13735, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 34708, 0, 3,
                                                                       31348, 12220, 31618, 3418,
                                                                       3473, 13900, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 35038, 0, 3,
                                                                       31618, 12355, 31888, 3473,
                                                                       3528, 14065, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 35368, 0, 3,
                                                                       31888, 12490, 32158, 3528,
                                                                       3583, 14230, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 35698, 0, 3,
                                                                       32428, 12760, 32698, 3693,
                                                                       3748, 14395, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 36028, 0, 3,
                                                                       32698, 12895, 32968, 3748,
                                                                       3803, 14560, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 36358, 0, 3,
                                                                       32968, 13030, 33238, 3803,
                                                                       3858, 14725, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 36688, 0, 3,
                                                                       33238, 13165, 33508, 3858,
                                                                       3913, 14890, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 37018, 0, 3,
                                                                       33508, 13300, 33778, 3913,
                                                                       3968, 15055, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 37348, 0, 3,
                                                                       34048, 13570, 34378, 4078,
                                                                       4144, 15220, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 37744, 0, 3,
                                                                       34378, 13735, 34708, 4144,
                                                                       4210, 15418, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 38140, 0, 3,
                                                                       34708, 13900, 35038, 4210,
                                                                       4276, 15616, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 38536, 0, 3,
                                                                       35038, 14065, 35368, 4276,
                                                                       4342, 15814, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 38932, 0, 3,
                                                                       35698, 14395, 36028, 4474,
                                                                       4540, 16012, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 39328, 0, 3,
                                                                       36028, 14560, 36358, 4540,
                                                                       4606, 16210, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 39724, 0, 3,
                                                                       36358, 14725, 36688, 4606,
                                                                       4672, 16408, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 40120, 0, 3,
                                                                       36688, 14890, 37018, 4672,
                                                                       4738, 16606, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 40516, 0, 3,
                                                                       37348, 15220, 37744, 4870,
                                                                       4948, 16804, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 40984, 0, 3,
                                                                       37744, 15418, 38140, 4948,
                                                                       5026, 17038, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 41452, 0, 3,
                                                                       38140, 15616, 38536, 5026,
                                                                       5104, 17272, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 41920, 0, 3,
                                                                       38932, 16012, 39328, 5260,
                                                                       5338, 17506, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 42388, 0, 3,
                                                                       39328, 16210, 39724, 5338,
                                                                       5416, 17740, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 42856, 0, 3,
                                                                       39724, 16408, 40120, 5416,
                                                                       5494, 17974, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43324, 3, 5650,
                                                                       5653, 18220, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43334, 3, 5653,
                                                                       5656, 18226, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43344, 3, 5656,
                                                                       5659, 18232, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43354, 3, 5659,
                                                                       5662, 18238, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43364, 3, 5662,
                                                                       5665, 18244, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43374, 3, 5665,
                                                                       5668, 18250, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43384, 3, 5668,
                                                                       5671, 18256, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43394, 3, 5671,
                                                                       5674, 18262, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43404, 3, 5674,
                                                                       5677, 18268, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43414, 3, 5677,
                                                                       5680, 18274, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43424, 3, 5680,
                                                                       5683, 18280, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43434, 3, 5683,
                                                                       5686, 18286, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43444, 3, 5692,
                                                                       5695, 18304, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43454, 3, 5695,
                                                                       5698, 18310, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43464, 3, 5698,
                                                                       5701, 18316, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43474, 3, 5701,
                                                                       5704, 18322, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43484, 3, 5704,
                                                                       5707, 18328, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43494, 3, 5707,
                                                                       5710, 18334, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43504, 3, 5710,
                                                                       5713, 18340, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43514, 3, 5713,
                                                                       5716, 18346, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43524, 3, 5716,
                                                                       5719, 18352, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43534, 3, 5719,
                                                                       5722, 18358, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43544, 3, 5722,
                                                                       5725, 18364, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43554, 3, 5725,
                                                                       5728, 18370, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 43564, 0, 3,
                                                                       43324, 18220, 43334, 5734,
                                                                       5743, 18412, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 43594, 0, 3,
                                                                       43334, 18226, 43344, 5743,
                                                                       5752, 18430, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 43624, 0, 3,
                                                                       43344, 18232, 43354, 5752,
                                                                       5761, 18448, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 43654, 0, 3,
                                                                       43354, 18238, 43364, 5761,
                                                                       5770, 18466, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 43684, 0, 3,
                                                                       43364, 18244, 43374, 5770,
                                                                       5779, 18484, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 43714, 0, 3,
                                                                       43374, 18250, 43384, 5779,
                                                                       5788, 18502, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 43744, 0, 3,
                                                                       43384, 18256, 43394, 5788,
                                                                       5797, 18520, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 43774, 0, 3,
                                                                       43394, 18262, 43404, 5797,
                                                                       5806, 18538, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 43804, 0, 3,
                                                                       43404, 18268, 43414, 5806,
                                                                       5815, 18556, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 43834, 0, 3,
                                                                       43414, 18274, 43424, 5815,
                                                                       5824, 18574, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 43864, 0, 3,
                                                                       43424, 18280, 43434, 5824,
                                                                       5833, 18592, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 43894, 0, 3,
                                                                       43444, 18304, 43454, 5851,
                                                                       5860, 18646, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 43924, 0, 3,
                                                                       43454, 18310, 43464, 5860,
                                                                       5869, 18664, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 43954, 0, 3,
                                                                       43464, 18316, 43474, 5869,
                                                                       5878, 18682, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 43984, 0, 3,
                                                                       43474, 18322, 43484, 5878,
                                                                       5887, 18700, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44014, 0, 3,
                                                                       43484, 18328, 43494, 5887,
                                                                       5896, 18718, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44044, 0, 3,
                                                                       43494, 18334, 43504, 5896,
                                                                       5905, 18736, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44074, 0, 3,
                                                                       43504, 18340, 43514, 5905,
                                                                       5914, 18754, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44104, 0, 3,
                                                                       43514, 18346, 43524, 5914,
                                                                       5923, 18772, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44134, 0, 3,
                                                                       43524, 18352, 43534, 5923,
                                                                       5932, 18790, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44164, 0, 3,
                                                                       43534, 18358, 43544, 5932,
                                                                       5941, 18808, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44194, 0, 3,
                                                                       43544, 18364, 43554, 5941,
                                                                       5950, 18826, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 44224, 0, 3,
                                                                       43564, 18412, 43594, 5968,
                                                                       5986, 18916, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 44284, 0, 3,
                                                                       43594, 18430, 43624, 5986,
                                                                       6004, 18952, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 44344, 0, 3,
                                                                       43624, 18448, 43654, 6004,
                                                                       6022, 18988, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 44404, 0, 3,
                                                                       43654, 18466, 43684, 6022,
                                                                       6040, 19024, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 44464, 0, 3,
                                                                       43684, 18484, 43714, 6040,
                                                                       6058, 19060, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 44524, 0, 3,
                                                                       43714, 18502, 43744, 6058,
                                                                       6076, 19096, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 44584, 0, 3,
                                                                       43744, 18520, 43774, 6076,
                                                                       6094, 19132, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 44644, 0, 3,
                                                                       43774, 18538, 43804, 6094,
                                                                       6112, 19168, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 44704, 0, 3,
                                                                       43804, 18556, 43834, 6112,
                                                                       6130, 19204, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 44764, 0, 3,
                                                                       43834, 18574, 43864, 6130,
                                                                       6148, 19240, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 44824, 0, 3,
                                                                       43894, 18646, 43924, 6184,
                                                                       6202, 19348, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 44884, 0, 3,
                                                                       43924, 18664, 43954, 6202,
                                                                       6220, 19384, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 44944, 0, 3,
                                                                       43954, 18682, 43984, 6220,
                                                                       6238, 19420, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45004, 0, 3,
                                                                       43984, 18700, 44014, 6238,
                                                                       6256, 19456, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45064, 0, 3,
                                                                       44014, 18718, 44044, 6256,
                                                                       6274, 19492, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45124, 0, 3,
                                                                       44044, 18736, 44074, 6274,
                                                                       6292, 19528, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45184, 0, 3,
                                                                       44074, 18754, 44104, 6292,
                                                                       6310, 19564, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45244, 0, 3,
                                                                       44104, 18772, 44134, 6310,
                                                                       6328, 19600, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45304, 0, 3,
                                                                       44134, 18790, 44164, 6328,
                                                                       6346, 19636, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45364, 0, 3,
                                                                       44164, 18808, 44194, 6346,
                                                                       6364, 19672, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 45424, 0, 3,
                                                                       44224, 18916, 44284, 6400,
                                                                       6430, 19828, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 45524, 0, 3,
                                                                       44284, 18952, 44344, 6430,
                                                                       6460, 19888, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 45624, 0, 3,
                                                                       44344, 18988, 44404, 6460,
                                                                       6490, 19948, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 45724, 0, 3,
                                                                       44404, 19024, 44464, 6490,
                                                                       6520, 20008, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 45824, 0, 3,
                                                                       44464, 19060, 44524, 6520,
                                                                       6550, 20068, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 45924, 0, 3,
                                                                       44524, 19096, 44584, 6550,
                                                                       6580, 20128, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46024, 0, 3,
                                                                       44584, 19132, 44644, 6580,
                                                                       6610, 20188, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46124, 0, 3,
                                                                       44644, 19168, 44704, 6610,
                                                                       6640, 20248, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46224, 0, 3,
                                                                       44704, 19204, 44764, 6640,
                                                                       6670, 20308, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46324, 0, 3,
                                                                       44824, 19348, 44884, 6730,
                                                                       6760, 20488, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46424, 0, 3,
                                                                       44884, 19384, 44944, 6760,
                                                                       6790, 20548, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46524, 0, 3,
                                                                       44944, 19420, 45004, 6790,
                                                                       6820, 20608, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46624, 0, 3,
                                                                       45004, 19456, 45064, 6820,
                                                                       6850, 20668, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46724, 0, 3,
                                                                       45064, 19492, 45124, 6850,
                                                                       6880, 20728, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46824, 0, 3,
                                                                       45124, 19528, 45184, 6880,
                                                                       6910, 20788, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46924, 0, 3,
                                                                       45184, 19564, 45244, 6910,
                                                                       6940, 20848, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47024, 0, 3,
                                                                       45244, 19600, 45304, 6940,
                                                                       6970, 20908, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47124, 0, 3,
                                                                       45304, 19636, 45364, 6970,
                                                                       7000, 20968, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 47224, 0, 3,
                                                                       45424, 19828, 45524, 7060,
                                                                       7105, 21208, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 47374, 0, 3,
                                                                       45524, 19888, 45624, 7105,
                                                                       7150, 21298, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 47524, 0, 3,
                                                                       45624, 19948, 45724, 7150,
                                                                       7195, 21388, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 47674, 0, 3,
                                                                       45724, 20008, 45824, 7195,
                                                                       7240, 21478, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 47824, 0, 3,
                                                                       45824, 20068, 45924, 7240,
                                                                       7285, 21568, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 47974, 0, 3,
                                                                       45924, 20128, 46024, 7285,
                                                                       7330, 21658, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48124, 0, 3,
                                                                       46024, 20188, 46124, 7330,
                                                                       7375, 21748, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48274, 0, 3,
                                                                       46124, 20248, 46224, 7375,
                                                                       7420, 21838, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48424, 0, 3,
                                                                       46324, 20488, 46424, 7510,
                                                                       7555, 22108, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48574, 0, 3,
                                                                       46424, 20548, 46524, 7555,
                                                                       7600, 22198, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48724, 0, 3,
                                                                       46524, 20608, 46624, 7600,
                                                                       7645, 22288, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48874, 0, 3,
                                                                       46624, 20668, 46724, 7645,
                                                                       7690, 22378, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49024, 0, 3,
                                                                       46724, 20728, 46824, 7690,
                                                                       7735, 22468, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49174, 0, 3,
                                                                       46824, 20788, 46924, 7735,
                                                                       7780, 22558, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49324, 0, 3,
                                                                       46924, 20848, 47024, 7780,
                                                                       7825, 22648, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49474, 0, 3,
                                                                       47024, 20908, 47124, 7825,
                                                                       7870, 22738, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 49624, 0, 3,
                                                                       47224, 21208, 47374, 7960,
                                                                       8023, 23080, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 49834, 0, 3,
                                                                       47374, 21298, 47524, 8023,
                                                                       8086, 23206, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 50044, 0, 3,
                                                                       47524, 21388, 47674, 8086,
                                                                       8149, 23332, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 50254, 0, 3,
                                                                       47674, 21478, 47824, 8149,
                                                                       8212, 23458, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 50464, 0, 3,
                                                                       47824, 21568, 47974, 8212,
                                                                       8275, 23584, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 50674, 0, 3,
                                                                       47974, 21658, 48124, 8275,
                                                                       8338, 23710, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 50884, 0, 3,
                                                                       48124, 21748, 48274, 8338,
                                                                       8401, 23836, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51094, 0, 3,
                                                                       48424, 22108, 48574, 8527,
                                                                       8590, 24214, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51304, 0, 3,
                                                                       48574, 22198, 48724, 8590,
                                                                       8653, 24340, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51514, 0, 3,
                                                                       48724, 22288, 48874, 8653,
                                                                       8716, 24466, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51724, 0, 3,
                                                                       48874, 22378, 49024, 8716,
                                                                       8779, 24592, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51934, 0, 3,
                                                                       49024, 22468, 49174, 8779,
                                                                       8842, 24718, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52144, 0, 3,
                                                                       49174, 22558, 49324, 8842,
                                                                       8905, 24844, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52354, 0, 3,
                                                                       49324, 22648, 49474, 8905,
                                                                       8968, 24970, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 52564, 0, 3,
                                                                       49624, 23080, 49834, 9094,
                                                                       9178, 25432, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 52844, 0, 3,
                                                                       49834, 23206, 50044, 9178,
                                                                       9262, 25600, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 53124, 0, 3,
                                                                       50044, 23332, 50254, 9262,
                                                                       9346, 25768, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 53404, 0, 3,
                                                                       50254, 23458, 50464, 9346,
                                                                       9430, 25936, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 53684, 0, 3,
                                                                       50464, 23584, 50674, 9430,
                                                                       9514, 26104, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 53964, 0, 3,
                                                                       50674, 23710, 50884, 9514,
                                                                       9598, 26272, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54244, 0, 3,
                                                                       51094, 24214, 51304, 9766,
                                                                       9850, 26776, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54524, 0, 3,
                                                                       51304, 24340, 51514, 9850,
                                                                       9934, 26944, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54804, 0, 3,
                                                                       51514, 24466, 51724, 9934,
                                                                       10018, 27112, ncols,
                                                                       gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 55084, 0, 3,
                                                                       51724, 24592, 51934,
                                                                       10018, 10102, 27280,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 55364, 0, 3,
                                                                       51934, 24718, 52144,
                                                                       10102, 10186, 27448,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 55644, 0, 3,
                                                                       52144, 24844, 52354,
                                                                       10186, 10270, 27616,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 55924, 0, 3,
                                                                       52564, 25432, 52844,
                                                                       10438, 10546, 28216,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 56284, 0, 3,
                                                                       52844, 25600, 53124,
                                                                       10546, 10654, 28432,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 56644, 0, 3,
                                                                       53124, 25768, 53404,
                                                                       10654, 10762, 28648,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 57004, 0, 3,
                                                                       53404, 25936, 53684,
                                                                       10762, 10870, 28864,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 57364, 0, 3,
                                                                       53684, 26104, 53964,
                                                                       10870, 10978, 29080,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 57724, 0, 3,
                                                                       54244, 26776, 54524,
                                                                       11194, 11302, 29728,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 58084, 0, 3,
                                                                       54524, 26944, 54804,
                                                                       11302, 11410, 29944,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 58444, 0, 3,
                                                                       54804, 27112, 55084,
                                                                       11410, 11518, 30160,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 58804, 0, 3,
                                                                       55084, 27280, 55364,
                                                                       11518, 11626, 30376,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 59164, 0, 3,
                                                                       55364, 27448, 55644,
                                                                       11626, 11734, 30592,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 59524, 0, 3,
                                                                       55924, 28216, 56284,
                                                                       11950, 12085, 31348,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 59974, 0, 3,
                                                                       56284, 28432, 56644,
                                                                       12085, 12220, 31618,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 60424, 0, 3,
                                                                       56644, 28648, 57004,
                                                                       12220, 12355, 31888,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 60874, 0, 3,
                                                                       57004, 28864, 57364,
                                                                       12355, 12490, 32158,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 61324, 0, 3,
                                                                       57724, 29728, 58084,
                                                                       12760, 12895, 32968,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 61774, 0, 3,
                                                                       58084, 29944, 58444,
                                                                       12895, 13030, 33238,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 62224, 0, 3,
                                                                       58444, 30160, 58804,
                                                                       13030, 13165, 33508,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 62674, 0, 3,
                                                                       58804, 30376, 59164,
                                                                       13165, 13300, 33778,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 63124, 0, 3,
                                                                       59524, 31348, 59974,
                                                                       13570, 13735, 34708,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 63674, 0, 3,
                                                                       59974, 31618, 60424,
                                                                       13735, 13900, 35038,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 64224, 0, 3,
                                                                       60424, 31888, 60874,
                                                                       13900, 14065, 35368,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 64774, 0, 3,
                                                                       61324, 32968, 61774,
                                                                       14395, 14560, 36358,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 65324, 0, 3,
                                                                       61774, 33238, 62224,
                                                                       14560, 14725, 36688,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 65874, 0, 3,
                                                                       62224, 33508, 62674,
                                                                       14725, 14890, 37018,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 66424, 0, 3,
                                                                       63124, 34708, 63674,
                                                                       15220, 15418, 38140,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 67084, 0, 3,
                                                                       63674, 35038, 64224,
                                                                       15418, 15616, 38536,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 67744, 0, 3,
                                                                       64774, 36358, 65324,
                                                                       16012, 16210, 39724,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 68404, 0, 3,
                                                                       65324, 36688, 65874,
                                                                       16210, 16408, 40120,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 69064, 0, 3,
                                                                       66424, 38140, 67084,
                                                                       16804, 17038, 41452,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 69844, 0, 3,
                                                                       67744, 39724, 68404,
                                                                       17506, 17740, 42856,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70624, 3, 18208,
                                                                       18214, 43324, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70639, 3, 18214,
                                                                       18220, 43334, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70654, 3, 18220,
                                                                       18226, 43344, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70669, 3, 18226,
                                                                       18232, 43354, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70684, 3, 18232,
                                                                       18238, 43364, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70699, 3, 18238,
                                                                       18244, 43374, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70714, 3, 18244,
                                                                       18250, 43384, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70729, 3, 18250,
                                                                       18256, 43394, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70744, 3, 18256,
                                                                       18262, 43404, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70759, 3, 18262,
                                                                       18268, 43414, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70774, 3, 18268,
                                                                       18274, 43424, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70789, 3, 18274,
                                                                       18280, 43434, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70804, 3, 18292,
                                                                       18298, 43444, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70819, 3, 18298,
                                                                       18304, 43454, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70834, 3, 18304,
                                                                       18310, 43464, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70849, 3, 18310,
                                                                       18316, 43474, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70864, 3, 18316,
                                                                       18322, 43484, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70879, 3, 18322,
                                                                       18328, 43494, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70894, 3, 18328,
                                                                       18334, 43504, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70909, 3, 18334,
                                                                       18340, 43514, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70924, 3, 18340,
                                                                       18346, 43524, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70939, 3, 18346,
                                                                       18352, 43534, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70954, 3, 18352,
                                                                       18358, 43544, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 70969, 3, 18358,
                                                                       18364, 43554, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 70984, 0, 3,
                                                                       70624, 43324, 70639,
                                                                       18376, 18394, 43564,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71029, 0, 3,
                                                                       70639, 43334, 70654,
                                                                       18394, 18412, 43594,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71074, 0, 3,
                                                                       70654, 43344, 70669,
                                                                       18412, 18430, 43624,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71119, 0, 3,
                                                                       70669, 43354, 70684,
                                                                       18430, 18448, 43654,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71164, 0, 3,
                                                                       70684, 43364, 70699,
                                                                       18448, 18466, 43684,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71209, 0, 3,
                                                                       70699, 43374, 70714,
                                                                       18466, 18484, 43714,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71254, 0, 3,
                                                                       70714, 43384, 70729,
                                                                       18484, 18502, 43744,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71299, 0, 3,
                                                                       70729, 43394, 70744,
                                                                       18502, 18520, 43774,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71344, 0, 3,
                                                                       70744, 43404, 70759,
                                                                       18520, 18538, 43804,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71389, 0, 3,
                                                                       70759, 43414, 70774,
                                                                       18538, 18556, 43834,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71434, 0, 3,
                                                                       70774, 43424, 70789,
                                                                       18556, 18574, 43864,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71479, 0, 3,
                                                                       70804, 43444, 70819,
                                                                       18610, 18628, 43894,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71524, 0, 3,
                                                                       70819, 43454, 70834,
                                                                       18628, 18646, 43924,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71569, 0, 3,
                                                                       70834, 43464, 70849,
                                                                       18646, 18664, 43954,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71614, 0, 3,
                                                                       70849, 43474, 70864,
                                                                       18664, 18682, 43984,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71659, 0, 3,
                                                                       70864, 43484, 70879,
                                                                       18682, 18700, 44014,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71704, 0, 3,
                                                                       70879, 43494, 70894,
                                                                       18700, 18718, 44044,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71749, 0, 3,
                                                                       70894, 43504, 70909,
                                                                       18718, 18736, 44074,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71794, 0, 3,
                                                                       70909, 43514, 70924,
                                                                       18736, 18754, 44104,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71839, 0, 3,
                                                                       70924, 43524, 70939,
                                                                       18754, 18772, 44134,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71884, 0, 3,
                                                                       70939, 43534, 70954,
                                                                       18772, 18790, 44164,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 71929, 0, 3,
                                                                       70954, 43544, 70969,
                                                                       18790, 18808, 44194,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 71974, 0, 3,
                                                                       70984, 43564, 71029,
                                                                       18844, 18880, 44224,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 72064, 0, 3,
                                                                       71029, 43594, 71074,
                                                                       18880, 18916, 44284,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 72154, 0, 3,
                                                                       71074, 43624, 71119,
                                                                       18916, 18952, 44344,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 72244, 0, 3,
                                                                       71119, 43654, 71164,
                                                                       18952, 18988, 44404,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 72334, 0, 3,
                                                                       71164, 43684, 71209,
                                                                       18988, 19024, 44464,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 72424, 0, 3,
                                                                       71209, 43714, 71254,
                                                                       19024, 19060, 44524,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 72514, 0, 3,
                                                                       71254, 43744, 71299,
                                                                       19060, 19096, 44584,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 72604, 0, 3,
                                                                       71299, 43774, 71344,
                                                                       19096, 19132, 44644,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 72694, 0, 3,
                                                                       71344, 43804, 71389,
                                                                       19132, 19168, 44704,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 72784, 0, 3,
                                                                       71389, 43834, 71434,
                                                                       19168, 19204, 44764,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 72874, 0, 3,
                                                                       71479, 43894, 71524,
                                                                       19276, 19312, 44824,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 72964, 0, 3,
                                                                       71524, 43924, 71569,
                                                                       19312, 19348, 44884,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 73054, 0, 3,
                                                                       71569, 43954, 71614,
                                                                       19348, 19384, 44944,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 73144, 0, 3,
                                                                       71614, 43984, 71659,
                                                                       19384, 19420, 45004,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 73234, 0, 3,
                                                                       71659, 44014, 71704,
                                                                       19420, 19456, 45064,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 73324, 0, 3,
                                                                       71704, 44044, 71749,
                                                                       19456, 19492, 45124,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 73414, 0, 3,
                                                                       71749, 44074, 71794,
                                                                       19492, 19528, 45184,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 73504, 0, 3,
                                                                       71794, 44104, 71839,
                                                                       19528, 19564, 45244,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 73594, 0, 3,
                                                                       71839, 44134, 71884,
                                                                       19564, 19600, 45304,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 73684, 0, 3,
                                                                       71884, 44164, 71929,
                                                                       19600, 19636, 45364,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 73774, 0, 3,
                                                                       71974, 44224, 72064,
                                                                       19708, 19768, 45424,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 73924, 0, 3,
                                                                       72064, 44284, 72154,
                                                                       19768, 19828, 45524,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 74074, 0, 3,
                                                                       72154, 44344, 72244,
                                                                       19828, 19888, 45624,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 74224, 0, 3,
                                                                       72244, 44404, 72334,
                                                                       19888, 19948, 45724,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 74374, 0, 3,
                                                                       72334, 44464, 72424,
                                                                       19948, 20008, 45824,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 74524, 0, 3,
                                                                       72424, 44524, 72514,
                                                                       20008, 20068, 45924,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 74674, 0, 3,
                                                                       72514, 44584, 72604,
                                                                       20068, 20128, 46024,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 74824, 0, 3,
                                                                       72604, 44644, 72694,
                                                                       20128, 20188, 46124,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 74974, 0, 3,
                                                                       72694, 44704, 72784,
                                                                       20188, 20248, 46224,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 75124, 0, 3,
                                                                       72874, 44824, 72964,
                                                                       20368, 20428, 46324,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 75274, 0, 3,
                                                                       72964, 44884, 73054,
                                                                       20428, 20488, 46424,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 75424, 0, 3,
                                                                       73054, 44944, 73144,
                                                                       20488, 20548, 46524,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 75574, 0, 3,
                                                                       73144, 45004, 73234,
                                                                       20548, 20608, 46624,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 75724, 0, 3,
                                                                       73234, 45064, 73324,
                                                                       20608, 20668, 46724,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 75874, 0, 3,
                                                                       73324, 45124, 73414,
                                                                       20668, 20728, 46824,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 76024, 0, 3,
                                                                       73414, 45184, 73504,
                                                                       20728, 20788, 46924,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 76174, 0, 3,
                                                                       73504, 45244, 73594,
                                                                       20788, 20848, 47024,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 76324, 0, 3,
                                                                       73594, 45304, 73684,
                                                                       20848, 20908, 47124,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 76474, 0, 3,
                                                                       73774, 45424, 73924,
                                                                       21028, 21118, 47224,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 76699, 0, 3,
                                                                       73924, 45524, 74074,
                                                                       21118, 21208, 47374,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 76924, 0, 3,
                                                                       74074, 45624, 74224,
                                                                       21208, 21298, 47524,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 77149, 0, 3,
                                                                       74224, 45724, 74374,
                                                                       21298, 21388, 47674,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 77374, 0, 3,
                                                                       74374, 45824, 74524,
                                                                       21388, 21478, 47824,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 77599, 0, 3,
                                                                       74524, 45924, 74674,
                                                                       21478, 21568, 47974,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 77824, 0, 3,
                                                                       74674, 46024, 74824,
                                                                       21568, 21658, 48124,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 78049, 0, 3,
                                                                       74824, 46124, 74974,
                                                                       21658, 21748, 48274,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 78274, 0, 3,
                                                                       75124, 46324, 75274,
                                                                       21928, 22018, 48424,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 78499, 0, 3,
                                                                       75274, 46424, 75424,
                                                                       22018, 22108, 48574,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 78724, 0, 3,
                                                                       75424, 46524, 75574,
                                                                       22108, 22198, 48724,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 78949, 0, 3,
                                                                       75574, 46624, 75724,
                                                                       22198, 22288, 48874,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 79174, 0, 3,
                                                                       75724, 46724, 75874,
                                                                       22288, 22378, 49024,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 79399, 0, 3,
                                                                       75874, 46824, 76024,
                                                                       22378, 22468, 49174,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 79624, 0, 3,
                                                                       76024, 46924, 76174,
                                                                       22468, 22558, 49324,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 79849, 0, 3,
                                                                       76174, 47024, 76324,
                                                                       22558, 22648, 49474,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 80074, 0, 3,
                                                                       76474, 47224, 76699,
                                                                       22828, 22954, 49624,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 80389, 0, 3,
                                                                       76699, 47374, 76924,
                                                                       22954, 23080, 49834,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 80704, 0, 3,
                                                                       76924, 47524, 77149,
                                                                       23080, 23206, 50044,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 81019, 0, 3,
                                                                       77149, 47674, 77374,
                                                                       23206, 23332, 50254,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 81334, 0, 3,
                                                                       77374, 47824, 77599,
                                                                       23332, 23458, 50464,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 81649, 0, 3,
                                                                       77599, 47974, 77824,
                                                                       23458, 23584, 50674,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 81964, 0, 3,
                                                                       77824, 48124, 78049,
                                                                       23584, 23710, 50884,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 82279, 0, 3,
                                                                       78274, 48424, 78499,
                                                                       23962, 24088, 51094,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 82594, 0, 3,
                                                                       78499, 48574, 78724,
                                                                       24088, 24214, 51304,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 82909, 0, 3,
                                                                       78724, 48724, 78949,
                                                                       24214, 24340, 51514,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 83224, 0, 3,
                                                                       78949, 48874, 79174,
                                                                       24340, 24466, 51724,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 83539, 0, 3,
                                                                       79174, 49024, 79399,
                                                                       24466, 24592, 51934,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 83854, 0, 3,
                                                                       79399, 49174, 79624,
                                                                       24592, 24718, 52144,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 84169, 0, 3,
                                                                       79624, 49324, 79849,
                                                                       24718, 24844, 52354,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 84484, 0, 3,
                                                                       80074, 49624, 80389,
                                                                       25096, 25264, 52564,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 84904, 0, 3,
                                                                       80389, 49834, 80704,
                                                                       25264, 25432, 52844,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 85324, 0, 3,
                                                                       80704, 50044, 81019,
                                                                       25432, 25600, 53124,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 85744, 0, 3,
                                                                       81019, 50254, 81334,
                                                                       25600, 25768, 53404,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 86164, 0, 3,
                                                                       81334, 50464, 81649,
                                                                       25768, 25936, 53684,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 86584, 0, 3,
                                                                       81649, 50674, 81964,
                                                                       25936, 26104, 53964,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 87004, 0, 3,
                                                                       82279, 51094, 82594,
                                                                       26440, 26608, 54244,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 87424, 0, 3,
                                                                       82594, 51304, 82909,
                                                                       26608, 26776, 54524,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 87844, 0, 3,
                                                                       82909, 51514, 83224,
                                                                       26776, 26944, 54804,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 88264, 0, 3,
                                                                       83224, 51724, 83539,
                                                                       26944, 27112, 55084,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 88684, 0, 3,
                                                                       83539, 51934, 83854,
                                                                       27112, 27280, 55364,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 89104, 0, 3,
                                                                       83854, 52144, 84169,
                                                                       27280, 27448, 55644,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 89524, 0, 3,
                                                                       84484, 52564, 84904,
                                                                       27784, 28000, 55924,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 90064, 0, 3,
                                                                       84904, 52844, 85324,
                                                                       28000, 28216, 56284,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 90604, 0, 3,
                                                                       85324, 53124, 85744,
                                                                       28216, 28432, 56644,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 91144, 0, 3,
                                                                       85744, 53404, 86164,
                                                                       28432, 28648, 57004,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 91684, 0, 3,
                                                                       86164, 53684, 86584,
                                                                       28648, 28864, 57364,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 92224, 0, 3,
                                                                       87004, 54244, 87424,
                                                                       29296, 29512, 57724,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 92764, 0, 3,
                                                                       87424, 54524, 87844,
                                                                       29512, 29728, 58084,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 93304, 0, 3,
                                                                       87844, 54804, 88264,
                                                                       29728, 29944, 58444,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 93844, 0, 3,
                                                                       88264, 55084, 88684,
                                                                       29944, 30160, 58804,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 94384, 0, 3,
                                                                       88684, 55364, 89104,
                                                                       30160, 30376, 59164,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 94924, 0, 3,
                                                                       89524, 55924, 90064,
                                                                       30808, 31078, 59524,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 95599, 0, 3,
                                                                       90064, 56284, 90604,
                                                                       31078, 31348, 59974,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 96274, 0, 3,
                                                                       90604, 56644, 91144,
                                                                       31348, 31618, 60424,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 96949, 0, 3,
                                                                       91144, 57004, 91684,
                                                                       31618, 31888, 60874,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 97624, 0, 3,
                                                                       92224, 57724, 92764,
                                                                       32428, 32698, 61324,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 98299, 0, 3,
                                                                       92764, 58084, 93304,
                                                                       32698, 32968, 61774,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 98974, 0, 3,
                                                                       93304, 58444, 93844,
                                                                       32968, 33238, 62224,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 99649, 0, 3,
                                                                       93844, 58804, 94384,
                                                                       33238, 33508, 62674,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 100324, 0, 3,
                                                                       94924, 59524, 95599,
                                                                       34048, 34378, 63124,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 101149, 0, 3,
                                                                       95599, 59974, 96274,
                                                                       34378, 34708, 63674,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 101974, 0, 3,
                                                                       96274, 60424, 96949,
                                                                       34708, 35038, 64224,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 102799, 0, 3,
                                                                       97624, 61324, 98299,
                                                                       35698, 36028, 64774,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 103624, 0, 3,
                                                                       98299, 61774, 98974,
                                                                       36028, 36358, 65324,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 104449, 0, 3,
                                                                       98974, 62224, 99649,
                                                                       36358, 36688, 65874,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 105274, 0, 3,
                                                                       100324, 63124, 101149,
                                                                       37348, 37744, 66424,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 106264, 0, 3,
                                                                       101149, 63674, 101974,
                                                                       37744, 38140, 67084,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 107254, 0, 3,
                                                                       102799, 64774, 103624,
                                                                       38932, 39328, 67744,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 108244, 0, 3,
                                                                       103624, 65324, 104449,
                                                                       39328, 39724, 68404,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 109234, 0, 3,
                                                                       105274, 66424, 106264,
                                                                       40516, 40984, 69064,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 110404, 0, 3,
                                                                       107254, 67744, 108244,
                                                                       41920, 42388, 69844,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 111574, 84484, 420, ncols);

                    simdfunc::contract_primitives(buffer, 112246, 87004, 420, ncols);

                    simdfunc::contract_primitives(buffer, 112918, 89524, 540, ncols);

                    simdfunc::contract_primitives(buffer, 113782, 92224, 540, ncols);

                    simdfunc::contract_primitives(buffer, 114646, 94924, 675, ncols);

                    simdfunc::contract_primitives(buffer, 115726, 97624, 675, ncols);

                    simdfunc::contract_primitives(buffer, 116806, 100324, 825, ncols);

                    simdfunc::contract_primitives(buffer, 118126, 102799, 825, ncols);

                    simdfunc::contract_primitives(buffer, 119446, 105274, 990, ncols);

                    simdfunc::contract_primitives(buffer, 121030, 107254, 990, ncols);

                    simdfunc::contract_primitives(buffer, 122614, 109234, 1170, ncols);

                    simdfunc::contract_primitives(buffer, 124486, 110404, 1170, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 111994, 111574, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 112666, 112246, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 113458, 112918, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 114322, 113782, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 115321, 114646, 45, 1, nmax);

        simdtrf::transform_g_inner(buffer, 116401, 115726, 45, 1, nmax);

        simdtrf::transform_g_inner(buffer, 117631, 116806, 55, 1, nmax);

        simdtrf::transform_g_inner(buffer, 118951, 118126, 55, 1, nmax);

        simdtrf::transform_g_inner(buffer, 120436, 119446, 66, 1, nmax);

        simdtrf::transform_g_inner(buffer, 122020, 121030, 66, 1, nmax);

        simdtrf::transform_g_inner(buffer, 123784, 122614, 78, 1, nmax);

        simdtrf::transform_g_inner(buffer, 125656, 124486, 78, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 126358, 111994, 113458, 9,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 127114, 112666, 114322, 9,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 127870, 113458, 115321, 9,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 128842, 114322, 116401, 9,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 129814, 115321, 117631, 9,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 131029, 116401, 118951, 9,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 132244, 117631, 120436, 9,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 133729, 118951, 122020, 9,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 135214, 120436, 123784, 9,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 136996, 122020, 125656, 9,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 138778, 126358, 127870, 9,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 140290, 127114, 128842, 9,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 141802, 127870, 129814, 9,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 143746, 128842, 131029, 9,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 145690, 129814, 132244, 9,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 148120, 131029, 133729, 9,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 150550, 132244, 135214, 9,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 153520, 133729, 136996, 9,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 156490, 138778, 141802, 9,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 159010, 140290, 143746, 9,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 161530, 141802, 145690, 9,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 164770, 143746, 148120, 9,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 168010, 145690, 150550, 9,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 172060, 148120, 153520, 9,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 176110, 156490, 161530, 9,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 179890, 159010, 164770, 9,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 183670, 161530, 168010, 9,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 188530, 164770, 172060, 9,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 193390, 176110, 183670, 9,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 198682, 179890, 188530, 9,
                                             nmax);

        simdtrf::transform_h_inner(buffer, 203974, 198682, 28, 9, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 203974, 99, nmax);

        simdtrf::transform_h_inner(buffer, 203974, 193390, 28, 9, nmax);

        simdtrf::transform_i_outer(values + 1287 * nvalues + n * npairs, nvalues, buffer, 203974,
                                   99, nmax);
    }

    for (size_t m = 0; m < 2574; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
