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


#include "SimdThreeCenterElectronRepulsionRsRecHGK.hpp"

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
#include "SimdTransferHD.hpp"
#include "SimdTransferHF.hpp"
#include "SimdTransferHG.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIF.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_hgk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_hgk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 314863, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2970 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 314863, 242008, 18045, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4078, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4081, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4084, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4087, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4090, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4093, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4096, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4099, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4102, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4105, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4108, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4111, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4114, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4117, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4120, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4123, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4126, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4129, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4132, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4135, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4138, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4141, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4144, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4147, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4150, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4153, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4156, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4159, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4162, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4165, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4168, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4171, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4174, 3, 7, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4183, 3, 8, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4192, 3, 9, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4201, 3, 10, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4210, 3, 11, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4219, 3, 12, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4228, 3, 13, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4237, 3, 14, 61,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4246, 3, 15, 64,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4255, 3, 16, 67,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4264, 3, 17, 70,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4273, 3, 18, 73,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4282, 3, 19, 76,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4291, 3, 20, 79,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4300, 3, 21, 82,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4309, 3, 24, 85,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4318, 3, 25, 88,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4327, 3, 26, 91,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4336, 3, 27, 94,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4345, 3, 28, 97,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4354, 3, 29, 100,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4363, 3, 30, 103,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4372, 3, 31, 106,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4381, 3, 32, 109,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4390, 3, 33, 112,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4399, 3, 34, 115,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4408, 3, 35, 118,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4417, 3, 36, 121,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4426, 3, 37, 124,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4435, 3, 38, 127,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4444, 3, 40, 130,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4462, 3, 43, 136,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4480, 3, 46, 142,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4498, 3, 49, 148,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4516, 3, 52, 154,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4534, 3, 55, 160,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4552, 3, 58, 166,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4570, 3, 61, 172,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4588, 3, 64, 178,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4606, 3, 67, 184,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4624, 3, 70, 190,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4642, 3, 73, 196,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4660, 3, 76, 202,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4678, 3, 79, 208,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4696, 3, 85, 214,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4714, 3, 88, 220,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4732, 3, 91, 226,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4750, 3, 94, 232,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4768, 3, 97, 238,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4786, 3, 100, 244,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4804, 3, 103, 250,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4822, 3, 106, 256,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4840, 3, 109, 262,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4858, 3, 112, 268,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4876, 3, 115, 274,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4894, 3, 118, 280,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4912, 3, 121, 286,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4930, 3, 124, 292,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4948, 3, 130, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4978, 3, 136, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5008, 3, 142, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5038, 3, 148, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5068, 3, 154, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5098, 3, 160, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5128, 3, 166, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5158, 3, 172, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5188, 3, 178, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5218, 3, 184, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5248, 3, 190, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5278, 3, 196, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5308, 3, 202, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5338, 3, 214, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5368, 3, 220, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5398, 3, 226, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5428, 3, 232, 458,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5458, 3, 238, 468,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5488, 3, 244, 478,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5518, 3, 250, 488,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5548, 3, 256, 498,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5578, 3, 262, 508,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5608, 3, 268, 518,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5638, 3, 274, 528,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5668, 3, 280, 538,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5698, 3, 286, 548,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5728, 3, 298, 558,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5773, 3, 308, 573,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5818, 3, 318, 588,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5863, 3, 328, 603,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5908, 3, 338, 618,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5953, 3, 348, 633,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5998, 3, 358, 648,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6043, 3, 368, 663,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6088, 3, 378, 678,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6133, 3, 388, 693,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6178, 3, 398, 708,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6223, 3, 408, 723,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6268, 3, 428, 738,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6313, 3, 438, 753,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6358, 3, 448, 768,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6403, 3, 458, 783,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6448, 3, 468, 798,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6493, 3, 478, 813,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6538, 3, 488, 828,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6583, 3, 498, 843,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6628, 3, 508, 858,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6673, 3, 518, 873,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6718, 3, 528, 888,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6763, 3, 538, 903,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6808, 3, 558, 918,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6871, 3, 573, 939,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6934, 3, 588, 960,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6997, 3, 603, 981,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7060, 3, 618,
                                                                       1002, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7123, 3, 633,
                                                                       1023, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7186, 3, 648,
                                                                       1044, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7249, 3, 663,
                                                                       1065, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7312, 3, 678,
                                                                       1086, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7375, 3, 693,
                                                                       1107, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7438, 3, 708,
                                                                       1128, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7501, 3, 738,
                                                                       1149, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7564, 3, 753,
                                                                       1170, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7627, 3, 768,
                                                                       1191, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7690, 3, 783,
                                                                       1212, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7753, 3, 798,
                                                                       1233, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7816, 3, 813,
                                                                       1254, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7879, 3, 828,
                                                                       1275, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7942, 3, 843,
                                                                       1296, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8005, 3, 858,
                                                                       1317, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8068, 3, 873,
                                                                       1338, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8131, 3, 888,
                                                                       1359, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8194, 3, 918,
                                                                       1380, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8278, 3, 939,
                                                                       1408, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8362, 3, 960,
                                                                       1436, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8446, 3, 981,
                                                                       1464, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8530, 3, 1002,
                                                                       1492, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8614, 3, 1023,
                                                                       1520, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8698, 3, 1044,
                                                                       1548, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8782, 3, 1065,
                                                                       1576, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8866, 3, 1086,
                                                                       1604, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8950, 3, 1107,
                                                                       1632, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9034, 3, 1149,
                                                                       1660, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9118, 3, 1170,
                                                                       1688, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9202, 3, 1191,
                                                                       1716, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9286, 3, 1212,
                                                                       1744, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9370, 3, 1233,
                                                                       1772, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9454, 3, 1254,
                                                                       1800, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9538, 3, 1275,
                                                                       1828, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9622, 3, 1296,
                                                                       1856, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9706, 3, 1317,
                                                                       1884, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9790, 3, 1338,
                                                                       1912, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9874, 3, 1380,
                                                                       1940, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9982, 3, 1408,
                                                                       1976, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10090, 3, 1436,
                                                                       2012, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10198, 3, 1464,
                                                                       2048, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10306, 3, 1492,
                                                                       2084, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10414, 3, 1520,
                                                                       2120, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10522, 3, 1548,
                                                                       2156, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10630, 3, 1576,
                                                                       2192, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10738, 3, 1604,
                                                                       2228, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10846, 3, 1660,
                                                                       2264, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10954, 3, 1688,
                                                                       2300, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11062, 3, 1716,
                                                                       2336, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11170, 3, 1744,
                                                                       2372, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11278, 3, 1772,
                                                                       2408, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11386, 3, 1800,
                                                                       2444, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11494, 3, 1828,
                                                                       2480, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11602, 3, 1856,
                                                                       2516, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11710, 3, 1884,
                                                                       2552, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11818, 3, 1940,
                                                                       2588, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11953, 3, 1976,
                                                                       2633, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12088, 3, 2012,
                                                                       2678, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12223, 3, 2048,
                                                                       2723, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12358, 3, 2084,
                                                                       2768, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12493, 3, 2120,
                                                                       2813, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12628, 3, 2156,
                                                                       2858, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12763, 3, 2192,
                                                                       2903, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12898, 3, 2264,
                                                                       2948, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13033, 3, 2300,
                                                                       2993, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13168, 3, 2336,
                                                                       3038, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13303, 3, 2372,
                                                                       3083, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13438, 3, 2408,
                                                                       3128, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13573, 3, 2444,
                                                                       3173, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13708, 3, 2480,
                                                                       3218, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13843, 3, 2516,
                                                                       3263, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13978, 3, 2588,
                                                                       3308, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14143, 3, 2633,
                                                                       3363, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14308, 3, 2678,
                                                                       3418, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14473, 3, 2723,
                                                                       3473, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14638, 3, 2768,
                                                                       3528, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14803, 3, 2813,
                                                                       3583, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14968, 3, 2858,
                                                                       3638, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15133, 3, 2948,
                                                                       3693, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15298, 3, 2993,
                                                                       3748, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15463, 3, 3038,
                                                                       3803, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15628, 3, 3083,
                                                                       3858, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15793, 3, 3128,
                                                                       3913, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15958, 3, 3173,
                                                                       3968, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16123, 3, 3218,
                                                                       4023, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16288, 3, 7, 8,
                                                                       4084, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16294, 3, 8, 9,
                                                                       4087, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16300, 3, 9, 10,
                                                                       4090, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16306, 3, 10, 11,
                                                                       4093, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16312, 3, 11, 12,
                                                                       4096, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16318, 3, 12, 13,
                                                                       4099, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16324, 3, 13, 14,
                                                                       4102, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16330, 3, 14, 15,
                                                                       4105, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16336, 3, 15, 16,
                                                                       4108, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16342, 3, 16, 17,
                                                                       4111, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16348, 3, 17, 18,
                                                                       4114, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16354, 3, 18, 19,
                                                                       4117, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16360, 3, 19, 20,
                                                                       4120, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16366, 3, 20, 21,
                                                                       4123, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16372, 3, 24, 25,
                                                                       4132, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16378, 3, 25, 26,
                                                                       4135, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16384, 3, 26, 27,
                                                                       4138, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16390, 3, 27, 28,
                                                                       4141, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16396, 3, 28, 29,
                                                                       4144, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16402, 3, 29, 30,
                                                                       4147, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16408, 3, 30, 31,
                                                                       4150, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16414, 3, 31, 32,
                                                                       4153, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16420, 3, 32, 33,
                                                                       4156, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16426, 3, 33, 34,
                                                                       4159, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16432, 3, 34, 35,
                                                                       4162, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16438, 3, 35, 36,
                                                                       4165, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16444, 3, 36, 37,
                                                                       4168, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16450, 3, 37, 38,
                                                                       4171, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16456, 0, 3,
                                                                       16288, 4084, 16294, 4192,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16474, 0, 3,
                                                                       16294, 4087, 16300, 4201,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16492, 0, 3,
                                                                       16300, 4090, 16306, 4210,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16510, 0, 3,
                                                                       16306, 4093, 16312, 4219,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16528, 0, 3,
                                                                       16312, 4096, 16318, 4228,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16546, 0, 3,
                                                                       16318, 4099, 16324, 4237,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16564, 0, 3,
                                                                       16324, 4102, 16330, 4246,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16582, 0, 3,
                                                                       16330, 4105, 16336, 4255,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16600, 0, 3,
                                                                       16336, 4108, 16342, 4264,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16618, 0, 3,
                                                                       16342, 4111, 16348, 4273,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16636, 0, 3,
                                                                       16348, 4114, 16354, 4282,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16654, 0, 3,
                                                                       16354, 4117, 16360, 4291,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16672, 0, 3,
                                                                       16360, 4120, 16366, 4300,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16690, 0, 3,
                                                                       16372, 4132, 16378, 4327,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16708, 0, 3,
                                                                       16378, 4135, 16384, 4336,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16726, 0, 3,
                                                                       16384, 4138, 16390, 4345,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16744, 0, 3,
                                                                       16390, 4141, 16396, 4354,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16762, 0, 3,
                                                                       16396, 4144, 16402, 4363,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16780, 0, 3,
                                                                       16402, 4147, 16408, 4372,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16798, 0, 3,
                                                                       16408, 4150, 16414, 4381,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16816, 0, 3,
                                                                       16414, 4153, 16420, 4390,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16834, 0, 3,
                                                                       16420, 4156, 16426, 4399,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16852, 0, 3,
                                                                       16426, 4159, 16432, 4408,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16870, 0, 3,
                                                                       16432, 4162, 16438, 4417,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16888, 0, 3,
                                                                       16438, 4165, 16444, 4426,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16906, 0, 3,
                                                                       16444, 4168, 16450, 4435,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16924, 0, 3,
                                                                       16456, 4192, 16474, 130,
                                                                       136, 4480, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16960, 0, 3,
                                                                       16474, 4201, 16492, 136,
                                                                       142, 4498, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16996, 0, 3,
                                                                       16492, 4210, 16510, 142,
                                                                       148, 4516, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17032, 0, 3,
                                                                       16510, 4219, 16528, 148,
                                                                       154, 4534, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17068, 0, 3,
                                                                       16528, 4228, 16546, 154,
                                                                       160, 4552, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17104, 0, 3,
                                                                       16546, 4237, 16564, 160,
                                                                       166, 4570, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17140, 0, 3,
                                                                       16564, 4246, 16582, 166,
                                                                       172, 4588, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17176, 0, 3,
                                                                       16582, 4255, 16600, 172,
                                                                       178, 4606, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17212, 0, 3,
                                                                       16600, 4264, 16618, 178,
                                                                       184, 4624, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17248, 0, 3,
                                                                       16618, 4273, 16636, 184,
                                                                       190, 4642, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17284, 0, 3,
                                                                       16636, 4282, 16654, 190,
                                                                       196, 4660, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17320, 0, 3,
                                                                       16654, 4291, 16672, 196,
                                                                       202, 4678, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17356, 0, 3,
                                                                       16690, 4327, 16708, 214,
                                                                       220, 4732, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17392, 0, 3,
                                                                       16708, 4336, 16726, 220,
                                                                       226, 4750, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17428, 0, 3,
                                                                       16726, 4345, 16744, 226,
                                                                       232, 4768, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17464, 0, 3,
                                                                       16744, 4354, 16762, 232,
                                                                       238, 4786, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17500, 0, 3,
                                                                       16762, 4363, 16780, 238,
                                                                       244, 4804, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17536, 0, 3,
                                                                       16780, 4372, 16798, 244,
                                                                       250, 4822, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17572, 0, 3,
                                                                       16798, 4381, 16816, 250,
                                                                       256, 4840, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17608, 0, 3,
                                                                       16816, 4390, 16834, 256,
                                                                       262, 4858, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17644, 0, 3,
                                                                       16834, 4399, 16852, 262,
                                                                       268, 4876, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17680, 0, 3,
                                                                       16852, 4408, 16870, 268,
                                                                       274, 4894, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17716, 0, 3,
                                                                       16870, 4417, 16888, 274,
                                                                       280, 4912, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17752, 0, 3,
                                                                       16888, 4426, 16906, 280,
                                                                       286, 4930, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 17788, 0, 3,
                                                                       16924, 4480, 16960, 298,
                                                                       308, 5008, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 17848, 0, 3,
                                                                       16960, 4498, 16996, 308,
                                                                       318, 5038, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 17908, 0, 3,
                                                                       16996, 4516, 17032, 318,
                                                                       328, 5068, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 17968, 0, 3,
                                                                       17032, 4534, 17068, 328,
                                                                       338, 5098, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18028, 0, 3,
                                                                       17068, 4552, 17104, 338,
                                                                       348, 5128, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18088, 0, 3,
                                                                       17104, 4570, 17140, 348,
                                                                       358, 5158, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18148, 0, 3,
                                                                       17140, 4588, 17176, 358,
                                                                       368, 5188, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18208, 0, 3,
                                                                       17176, 4606, 17212, 368,
                                                                       378, 5218, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18268, 0, 3,
                                                                       17212, 4624, 17248, 378,
                                                                       388, 5248, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18328, 0, 3,
                                                                       17248, 4642, 17284, 388,
                                                                       398, 5278, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18388, 0, 3,
                                                                       17284, 4660, 17320, 398,
                                                                       408, 5308, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18448, 0, 3,
                                                                       17356, 4732, 17392, 428,
                                                                       438, 5398, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18508, 0, 3,
                                                                       17392, 4750, 17428, 438,
                                                                       448, 5428, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18568, 0, 3,
                                                                       17428, 4768, 17464, 448,
                                                                       458, 5458, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18628, 0, 3,
                                                                       17464, 4786, 17500, 458,
                                                                       468, 5488, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18688, 0, 3,
                                                                       17500, 4804, 17536, 468,
                                                                       478, 5518, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18748, 0, 3,
                                                                       17536, 4822, 17572, 478,
                                                                       488, 5548, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18808, 0, 3,
                                                                       17572, 4840, 17608, 488,
                                                                       498, 5578, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18868, 0, 3,
                                                                       17608, 4858, 17644, 498,
                                                                       508, 5608, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18928, 0, 3,
                                                                       17644, 4876, 17680, 508,
                                                                       518, 5638, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18988, 0, 3,
                                                                       17680, 4894, 17716, 518,
                                                                       528, 5668, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19048, 0, 3,
                                                                       17716, 4912, 17752, 528,
                                                                       538, 5698, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19108, 0, 3,
                                                                       17788, 5008, 17848, 558,
                                                                       573, 5818, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19198, 0, 3,
                                                                       17848, 5038, 17908, 573,
                                                                       588, 5863, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19288, 0, 3,
                                                                       17908, 5068, 17968, 588,
                                                                       603, 5908, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19378, 0, 3,
                                                                       17968, 5098, 18028, 603,
                                                                       618, 5953, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19468, 0, 3,
                                                                       18028, 5128, 18088, 618,
                                                                       633, 5998, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19558, 0, 3,
                                                                       18088, 5158, 18148, 633,
                                                                       648, 6043, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19648, 0, 3,
                                                                       18148, 5188, 18208, 648,
                                                                       663, 6088, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19738, 0, 3,
                                                                       18208, 5218, 18268, 663,
                                                                       678, 6133, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19828, 0, 3,
                                                                       18268, 5248, 18328, 678,
                                                                       693, 6178, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19918, 0, 3,
                                                                       18328, 5278, 18388, 693,
                                                                       708, 6223, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20008, 0, 3,
                                                                       18448, 5398, 18508, 738,
                                                                       753, 6358, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20098, 0, 3,
                                                                       18508, 5428, 18568, 753,
                                                                       768, 6403, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20188, 0, 3,
                                                                       18568, 5458, 18628, 768,
                                                                       783, 6448, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20278, 0, 3,
                                                                       18628, 5488, 18688, 783,
                                                                       798, 6493, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20368, 0, 3,
                                                                       18688, 5518, 18748, 798,
                                                                       813, 6538, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20458, 0, 3,
                                                                       18748, 5548, 18808, 813,
                                                                       828, 6583, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20548, 0, 3,
                                                                       18808, 5578, 18868, 828,
                                                                       843, 6628, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20638, 0, 3,
                                                                       18868, 5608, 18928, 843,
                                                                       858, 6673, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20728, 0, 3,
                                                                       18928, 5638, 18988, 858,
                                                                       873, 6718, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20818, 0, 3,
                                                                       18988, 5668, 19048, 873,
                                                                       888, 6763, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 20908, 0, 3,
                                                                       19108, 5818, 19198, 918,
                                                                       939, 6934, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21034, 0, 3,
                                                                       19198, 5863, 19288, 939,
                                                                       960, 6997, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21160, 0, 3,
                                                                       19288, 5908, 19378, 960,
                                                                       981, 7060, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21286, 0, 3,
                                                                       19378, 5953, 19468, 981,
                                                                       1002, 7123, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21412, 0, 3,
                                                                       19468, 5998, 19558, 1002,
                                                                       1023, 7186, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21538, 0, 3,
                                                                       19558, 6043, 19648, 1023,
                                                                       1044, 7249, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21664, 0, 3,
                                                                       19648, 6088, 19738, 1044,
                                                                       1065, 7312, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21790, 0, 3,
                                                                       19738, 6133, 19828, 1065,
                                                                       1086, 7375, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21916, 0, 3,
                                                                       19828, 6178, 19918, 1086,
                                                                       1107, 7438, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22042, 0, 3,
                                                                       20008, 6358, 20098, 1149,
                                                                       1170, 7627, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22168, 0, 3,
                                                                       20098, 6403, 20188, 1170,
                                                                       1191, 7690, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22294, 0, 3,
                                                                       20188, 6448, 20278, 1191,
                                                                       1212, 7753, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22420, 0, 3,
                                                                       20278, 6493, 20368, 1212,
                                                                       1233, 7816, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22546, 0, 3,
                                                                       20368, 6538, 20458, 1233,
                                                                       1254, 7879, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22672, 0, 3,
                                                                       20458, 6583, 20548, 1254,
                                                                       1275, 7942, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22798, 0, 3,
                                                                       20548, 6628, 20638, 1275,
                                                                       1296, 8005, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22924, 0, 3,
                                                                       20638, 6673, 20728, 1296,
                                                                       1317, 8068, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23050, 0, 3,
                                                                       20728, 6718, 20818, 1317,
                                                                       1338, 8131, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23176, 0, 3,
                                                                       20908, 6934, 21034, 1380,
                                                                       1408, 8362, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23344, 0, 3,
                                                                       21034, 6997, 21160, 1408,
                                                                       1436, 8446, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23512, 0, 3,
                                                                       21160, 7060, 21286, 1436,
                                                                       1464, 8530, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23680, 0, 3,
                                                                       21286, 7123, 21412, 1464,
                                                                       1492, 8614, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23848, 0, 3,
                                                                       21412, 7186, 21538, 1492,
                                                                       1520, 8698, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24016, 0, 3,
                                                                       21538, 7249, 21664, 1520,
                                                                       1548, 8782, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24184, 0, 3,
                                                                       21664, 7312, 21790, 1548,
                                                                       1576, 8866, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24352, 0, 3,
                                                                       21790, 7375, 21916, 1576,
                                                                       1604, 8950, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24520, 0, 3,
                                                                       22042, 7627, 22168, 1660,
                                                                       1688, 9202, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24688, 0, 3,
                                                                       22168, 7690, 22294, 1688,
                                                                       1716, 9286, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24856, 0, 3,
                                                                       22294, 7753, 22420, 1716,
                                                                       1744, 9370, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25024, 0, 3,
                                                                       22420, 7816, 22546, 1744,
                                                                       1772, 9454, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25192, 0, 3,
                                                                       22546, 7879, 22672, 1772,
                                                                       1800, 9538, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25360, 0, 3,
                                                                       22672, 7942, 22798, 1800,
                                                                       1828, 9622, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25528, 0, 3,
                                                                       22798, 8005, 22924, 1828,
                                                                       1856, 9706, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25696, 0, 3,
                                                                       22924, 8068, 23050, 1856,
                                                                       1884, 9790, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25864, 0, 3,
                                                                       23176, 8362, 23344, 1940,
                                                                       1976, 10090, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26080, 0, 3,
                                                                       23344, 8446, 23512, 1976,
                                                                       2012, 10198, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26296, 0, 3,
                                                                       23512, 8530, 23680, 2012,
                                                                       2048, 10306, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26512, 0, 3,
                                                                       23680, 8614, 23848, 2048,
                                                                       2084, 10414, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26728, 0, 3,
                                                                       23848, 8698, 24016, 2084,
                                                                       2120, 10522, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26944, 0, 3,
                                                                       24016, 8782, 24184, 2120,
                                                                       2156, 10630, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27160, 0, 3,
                                                                       24184, 8866, 24352, 2156,
                                                                       2192, 10738, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27376, 0, 3,
                                                                       24520, 9202, 24688, 2264,
                                                                       2300, 11062, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27592, 0, 3,
                                                                       24688, 9286, 24856, 2300,
                                                                       2336, 11170, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27808, 0, 3,
                                                                       24856, 9370, 25024, 2336,
                                                                       2372, 11278, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28024, 0, 3,
                                                                       25024, 9454, 25192, 2372,
                                                                       2408, 11386, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28240, 0, 3,
                                                                       25192, 9538, 25360, 2408,
                                                                       2444, 11494, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28456, 0, 3,
                                                                       25360, 9622, 25528, 2444,
                                                                       2480, 11602, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28672, 0, 3,
                                                                       25528, 9706, 25696, 2480,
                                                                       2516, 11710, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28888, 0, 3,
                                                                       25864, 10090, 26080, 2588,
                                                                       2633, 12088, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29158, 0, 3,
                                                                       26080, 10198, 26296, 2633,
                                                                       2678, 12223, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29428, 0, 3,
                                                                       26296, 10306, 26512, 2678,
                                                                       2723, 12358, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29698, 0, 3,
                                                                       26512, 10414, 26728, 2723,
                                                                       2768, 12493, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29968, 0, 3,
                                                                       26728, 10522, 26944, 2768,
                                                                       2813, 12628, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 30238, 0, 3,
                                                                       26944, 10630, 27160, 2813,
                                                                       2858, 12763, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 30508, 0, 3,
                                                                       27376, 11062, 27592, 2948,
                                                                       2993, 13168, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 30778, 0, 3,
                                                                       27592, 11170, 27808, 2993,
                                                                       3038, 13303, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 31048, 0, 3,
                                                                       27808, 11278, 28024, 3038,
                                                                       3083, 13438, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 31318, 0, 3,
                                                                       28024, 11386, 28240, 3083,
                                                                       3128, 13573, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 31588, 0, 3,
                                                                       28240, 11494, 28456, 3128,
                                                                       3173, 13708, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 31858, 0, 3,
                                                                       28456, 11602, 28672, 3173,
                                                                       3218, 13843, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 32128, 0, 3,
                                                                       28888, 12088, 29158, 3308,
                                                                       3363, 14308, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 32458, 0, 3,
                                                                       29158, 12223, 29428, 3363,
                                                                       3418, 14473, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 32788, 0, 3,
                                                                       29428, 12358, 29698, 3418,
                                                                       3473, 14638, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 33118, 0, 3,
                                                                       29698, 12493, 29968, 3473,
                                                                       3528, 14803, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 33448, 0, 3,
                                                                       29968, 12628, 30238, 3528,
                                                                       3583, 14968, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 33778, 0, 3,
                                                                       30508, 13168, 30778, 3693,
                                                                       3748, 15463, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 34108, 0, 3,
                                                                       30778, 13303, 31048, 3748,
                                                                       3803, 15628, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 34438, 0, 3,
                                                                       31048, 13438, 31318, 3803,
                                                                       3858, 15793, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 34768, 0, 3,
                                                                       31318, 13573, 31588, 3858,
                                                                       3913, 15958, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 35098, 0, 3,
                                                                       31588, 13708, 31858, 3913,
                                                                       3968, 16123, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35428, 3, 4078,
                                                                       4081, 16288, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35438, 3, 4081,
                                                                       4084, 16294, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35448, 3, 4084,
                                                                       4087, 16300, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35458, 3, 4087,
                                                                       4090, 16306, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35468, 3, 4090,
                                                                       4093, 16312, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35478, 3, 4093,
                                                                       4096, 16318, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35488, 3, 4096,
                                                                       4099, 16324, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35498, 3, 4099,
                                                                       4102, 16330, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35508, 3, 4102,
                                                                       4105, 16336, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35518, 3, 4105,
                                                                       4108, 16342, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35528, 3, 4108,
                                                                       4111, 16348, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35538, 3, 4111,
                                                                       4114, 16354, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35548, 3, 4114,
                                                                       4117, 16360, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35558, 3, 4117,
                                                                       4120, 16366, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35568, 3, 4126,
                                                                       4129, 16372, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35578, 3, 4129,
                                                                       4132, 16378, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35588, 3, 4132,
                                                                       4135, 16384, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35598, 3, 4135,
                                                                       4138, 16390, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35608, 3, 4138,
                                                                       4141, 16396, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35618, 3, 4141,
                                                                       4144, 16402, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35628, 3, 4144,
                                                                       4147, 16408, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35638, 3, 4147,
                                                                       4150, 16414, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35648, 3, 4150,
                                                                       4153, 16420, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35658, 3, 4153,
                                                                       4156, 16426, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35668, 3, 4156,
                                                                       4159, 16432, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35678, 3, 4159,
                                                                       4162, 16438, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35688, 3, 4162,
                                                                       4165, 16444, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35698, 3, 4165,
                                                                       4168, 16450, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 35708, 0, 3,
                                                                       35428, 16288, 35438, 4174,
                                                                       4183, 16456, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 35738, 0, 3,
                                                                       35438, 16294, 35448, 4183,
                                                                       4192, 16474, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 35768, 0, 3,
                                                                       35448, 16300, 35458, 4192,
                                                                       4201, 16492, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 35798, 0, 3,
                                                                       35458, 16306, 35468, 4201,
                                                                       4210, 16510, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 35828, 0, 3,
                                                                       35468, 16312, 35478, 4210,
                                                                       4219, 16528, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 35858, 0, 3,
                                                                       35478, 16318, 35488, 4219,
                                                                       4228, 16546, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 35888, 0, 3,
                                                                       35488, 16324, 35498, 4228,
                                                                       4237, 16564, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 35918, 0, 3,
                                                                       35498, 16330, 35508, 4237,
                                                                       4246, 16582, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 35948, 0, 3,
                                                                       35508, 16336, 35518, 4246,
                                                                       4255, 16600, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 35978, 0, 3,
                                                                       35518, 16342, 35528, 4255,
                                                                       4264, 16618, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36008, 0, 3,
                                                                       35528, 16348, 35538, 4264,
                                                                       4273, 16636, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36038, 0, 3,
                                                                       35538, 16354, 35548, 4273,
                                                                       4282, 16654, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36068, 0, 3,
                                                                       35548, 16360, 35558, 4282,
                                                                       4291, 16672, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36098, 0, 3,
                                                                       35568, 16372, 35578, 4309,
                                                                       4318, 16690, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36128, 0, 3,
                                                                       35578, 16378, 35588, 4318,
                                                                       4327, 16708, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36158, 0, 3,
                                                                       35588, 16384, 35598, 4327,
                                                                       4336, 16726, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36188, 0, 3,
                                                                       35598, 16390, 35608, 4336,
                                                                       4345, 16744, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36218, 0, 3,
                                                                       35608, 16396, 35618, 4345,
                                                                       4354, 16762, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36248, 0, 3,
                                                                       35618, 16402, 35628, 4354,
                                                                       4363, 16780, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36278, 0, 3,
                                                                       35628, 16408, 35638, 4363,
                                                                       4372, 16798, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36308, 0, 3,
                                                                       35638, 16414, 35648, 4372,
                                                                       4381, 16816, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36338, 0, 3,
                                                                       35648, 16420, 35658, 4381,
                                                                       4390, 16834, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36368, 0, 3,
                                                                       35658, 16426, 35668, 4390,
                                                                       4399, 16852, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36398, 0, 3,
                                                                       35668, 16432, 35678, 4399,
                                                                       4408, 16870, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36428, 0, 3,
                                                                       35678, 16438, 35688, 4408,
                                                                       4417, 16888, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36458, 0, 3,
                                                                       35688, 16444, 35698, 4417,
                                                                       4426, 16906, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 36488, 0, 3,
                                                                       35708, 16456, 35738, 4444,
                                                                       4462, 16924, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 36548, 0, 3,
                                                                       35738, 16474, 35768, 4462,
                                                                       4480, 16960, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 36608, 0, 3,
                                                                       35768, 16492, 35798, 4480,
                                                                       4498, 16996, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 36668, 0, 3,
                                                                       35798, 16510, 35828, 4498,
                                                                       4516, 17032, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 36728, 0, 3,
                                                                       35828, 16528, 35858, 4516,
                                                                       4534, 17068, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 36788, 0, 3,
                                                                       35858, 16546, 35888, 4534,
                                                                       4552, 17104, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 36848, 0, 3,
                                                                       35888, 16564, 35918, 4552,
                                                                       4570, 17140, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 36908, 0, 3,
                                                                       35918, 16582, 35948, 4570,
                                                                       4588, 17176, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 36968, 0, 3,
                                                                       35948, 16600, 35978, 4588,
                                                                       4606, 17212, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37028, 0, 3,
                                                                       35978, 16618, 36008, 4606,
                                                                       4624, 17248, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37088, 0, 3,
                                                                       36008, 16636, 36038, 4624,
                                                                       4642, 17284, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37148, 0, 3,
                                                                       36038, 16654, 36068, 4642,
                                                                       4660, 17320, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37208, 0, 3,
                                                                       36098, 16690, 36128, 4696,
                                                                       4714, 17356, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37268, 0, 3,
                                                                       36128, 16708, 36158, 4714,
                                                                       4732, 17392, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37328, 0, 3,
                                                                       36158, 16726, 36188, 4732,
                                                                       4750, 17428, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37388, 0, 3,
                                                                       36188, 16744, 36218, 4750,
                                                                       4768, 17464, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37448, 0, 3,
                                                                       36218, 16762, 36248, 4768,
                                                                       4786, 17500, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37508, 0, 3,
                                                                       36248, 16780, 36278, 4786,
                                                                       4804, 17536, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37568, 0, 3,
                                                                       36278, 16798, 36308, 4804,
                                                                       4822, 17572, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37628, 0, 3,
                                                                       36308, 16816, 36338, 4822,
                                                                       4840, 17608, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37688, 0, 3,
                                                                       36338, 16834, 36368, 4840,
                                                                       4858, 17644, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37748, 0, 3,
                                                                       36368, 16852, 36398, 4858,
                                                                       4876, 17680, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37808, 0, 3,
                                                                       36398, 16870, 36428, 4876,
                                                                       4894, 17716, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37868, 0, 3,
                                                                       36428, 16888, 36458, 4894,
                                                                       4912, 17752, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 37928, 0, 3,
                                                                       36488, 16924, 36548, 4948,
                                                                       4978, 17788, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38028, 0, 3,
                                                                       36548, 16960, 36608, 4978,
                                                                       5008, 17848, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38128, 0, 3,
                                                                       36608, 16996, 36668, 5008,
                                                                       5038, 17908, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38228, 0, 3,
                                                                       36668, 17032, 36728, 5038,
                                                                       5068, 17968, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38328, 0, 3,
                                                                       36728, 17068, 36788, 5068,
                                                                       5098, 18028, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38428, 0, 3,
                                                                       36788, 17104, 36848, 5098,
                                                                       5128, 18088, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38528, 0, 3,
                                                                       36848, 17140, 36908, 5128,
                                                                       5158, 18148, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38628, 0, 3,
                                                                       36908, 17176, 36968, 5158,
                                                                       5188, 18208, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38728, 0, 3,
                                                                       36968, 17212, 37028, 5188,
                                                                       5218, 18268, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38828, 0, 3,
                                                                       37028, 17248, 37088, 5218,
                                                                       5248, 18328, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38928, 0, 3,
                                                                       37088, 17284, 37148, 5248,
                                                                       5278, 18388, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39028, 0, 3,
                                                                       37208, 17356, 37268, 5338,
                                                                       5368, 18448, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39128, 0, 3,
                                                                       37268, 17392, 37328, 5368,
                                                                       5398, 18508, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39228, 0, 3,
                                                                       37328, 17428, 37388, 5398,
                                                                       5428, 18568, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39328, 0, 3,
                                                                       37388, 17464, 37448, 5428,
                                                                       5458, 18628, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39428, 0, 3,
                                                                       37448, 17500, 37508, 5458,
                                                                       5488, 18688, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39528, 0, 3,
                                                                       37508, 17536, 37568, 5488,
                                                                       5518, 18748, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39628, 0, 3,
                                                                       37568, 17572, 37628, 5518,
                                                                       5548, 18808, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39728, 0, 3,
                                                                       37628, 17608, 37688, 5548,
                                                                       5578, 18868, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39828, 0, 3,
                                                                       37688, 17644, 37748, 5578,
                                                                       5608, 18928, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39928, 0, 3,
                                                                       37748, 17680, 37808, 5608,
                                                                       5638, 18988, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 40028, 0, 3,
                                                                       37808, 17716, 37868, 5638,
                                                                       5668, 19048, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40128, 0, 3,
                                                                       37928, 17788, 38028, 5728,
                                                                       5773, 19108, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40278, 0, 3,
                                                                       38028, 17848, 38128, 5773,
                                                                       5818, 19198, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40428, 0, 3,
                                                                       38128, 17908, 38228, 5818,
                                                                       5863, 19288, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40578, 0, 3,
                                                                       38228, 17968, 38328, 5863,
                                                                       5908, 19378, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40728, 0, 3,
                                                                       38328, 18028, 38428, 5908,
                                                                       5953, 19468, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40878, 0, 3,
                                                                       38428, 18088, 38528, 5953,
                                                                       5998, 19558, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 41028, 0, 3,
                                                                       38528, 18148, 38628, 5998,
                                                                       6043, 19648, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 41178, 0, 3,
                                                                       38628, 18208, 38728, 6043,
                                                                       6088, 19738, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 41328, 0, 3,
                                                                       38728, 18268, 38828, 6088,
                                                                       6133, 19828, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 41478, 0, 3,
                                                                       38828, 18328, 38928, 6133,
                                                                       6178, 19918, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 41628, 0, 3,
                                                                       39028, 18448, 39128, 6268,
                                                                       6313, 20008, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 41778, 0, 3,
                                                                       39128, 18508, 39228, 6313,
                                                                       6358, 20098, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 41928, 0, 3,
                                                                       39228, 18568, 39328, 6358,
                                                                       6403, 20188, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 42078, 0, 3,
                                                                       39328, 18628, 39428, 6403,
                                                                       6448, 20278, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 42228, 0, 3,
                                                                       39428, 18688, 39528, 6448,
                                                                       6493, 20368, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 42378, 0, 3,
                                                                       39528, 18748, 39628, 6493,
                                                                       6538, 20458, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 42528, 0, 3,
                                                                       39628, 18808, 39728, 6538,
                                                                       6583, 20548, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 42678, 0, 3,
                                                                       39728, 18868, 39828, 6583,
                                                                       6628, 20638, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 42828, 0, 3,
                                                                       39828, 18928, 39928, 6628,
                                                                       6673, 20728, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 42978, 0, 3,
                                                                       39928, 18988, 40028, 6673,
                                                                       6718, 20818, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 43128, 0, 3,
                                                                       40128, 19108, 40278, 6808,
                                                                       6871, 20908, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 43338, 0, 3,
                                                                       40278, 19198, 40428, 6871,
                                                                       6934, 21034, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 43548, 0, 3,
                                                                       40428, 19288, 40578, 6934,
                                                                       6997, 21160, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 43758, 0, 3,
                                                                       40578, 19378, 40728, 6997,
                                                                       7060, 21286, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 43968, 0, 3,
                                                                       40728, 19468, 40878, 7060,
                                                                       7123, 21412, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 44178, 0, 3,
                                                                       40878, 19558, 41028, 7123,
                                                                       7186, 21538, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 44388, 0, 3,
                                                                       41028, 19648, 41178, 7186,
                                                                       7249, 21664, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 44598, 0, 3,
                                                                       41178, 19738, 41328, 7249,
                                                                       7312, 21790, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 44808, 0, 3,
                                                                       41328, 19828, 41478, 7312,
                                                                       7375, 21916, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 45018, 0, 3,
                                                                       41628, 20008, 41778, 7501,
                                                                       7564, 22042, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 45228, 0, 3,
                                                                       41778, 20098, 41928, 7564,
                                                                       7627, 22168, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 45438, 0, 3,
                                                                       41928, 20188, 42078, 7627,
                                                                       7690, 22294, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 45648, 0, 3,
                                                                       42078, 20278, 42228, 7690,
                                                                       7753, 22420, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 45858, 0, 3,
                                                                       42228, 20368, 42378, 7753,
                                                                       7816, 22546, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 46068, 0, 3,
                                                                       42378, 20458, 42528, 7816,
                                                                       7879, 22672, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 46278, 0, 3,
                                                                       42528, 20548, 42678, 7879,
                                                                       7942, 22798, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 46488, 0, 3,
                                                                       42678, 20638, 42828, 7942,
                                                                       8005, 22924, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 46698, 0, 3,
                                                                       42828, 20728, 42978, 8005,
                                                                       8068, 23050, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 46908, 0, 3,
                                                                       43128, 20908, 43338, 8194,
                                                                       8278, 23176, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 47188, 0, 3,
                                                                       43338, 21034, 43548, 8278,
                                                                       8362, 23344, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 47468, 0, 3,
                                                                       43548, 21160, 43758, 8362,
                                                                       8446, 23512, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 47748, 0, 3,
                                                                       43758, 21286, 43968, 8446,
                                                                       8530, 23680, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 48028, 0, 3,
                                                                       43968, 21412, 44178, 8530,
                                                                       8614, 23848, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 48308, 0, 3,
                                                                       44178, 21538, 44388, 8614,
                                                                       8698, 24016, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 48588, 0, 3,
                                                                       44388, 21664, 44598, 8698,
                                                                       8782, 24184, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 48868, 0, 3,
                                                                       44598, 21790, 44808, 8782,
                                                                       8866, 24352, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 49148, 0, 3,
                                                                       45018, 22042, 45228, 9034,
                                                                       9118, 24520, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 49428, 0, 3,
                                                                       45228, 22168, 45438, 9118,
                                                                       9202, 24688, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 49708, 0, 3,
                                                                       45438, 22294, 45648, 9202,
                                                                       9286, 24856, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 49988, 0, 3,
                                                                       45648, 22420, 45858, 9286,
                                                                       9370, 25024, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 50268, 0, 3,
                                                                       45858, 22546, 46068, 9370,
                                                                       9454, 25192, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 50548, 0, 3,
                                                                       46068, 22672, 46278, 9454,
                                                                       9538, 25360, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 50828, 0, 3,
                                                                       46278, 22798, 46488, 9538,
                                                                       9622, 25528, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 51108, 0, 3,
                                                                       46488, 22924, 46698, 9622,
                                                                       9706, 25696, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 51388, 0, 3,
                                                                       46908, 23176, 47188, 9874,
                                                                       9982, 25864, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 51748, 0, 3,
                                                                       47188, 23344, 47468, 9982,
                                                                       10090, 26080, ncols,
                                                                       gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 52108, 0, 3,
                                                                       47468, 23512, 47748,
                                                                       10090, 10198, 26296,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 52468, 0, 3,
                                                                       47748, 23680, 48028,
                                                                       10198, 10306, 26512,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 52828, 0, 3,
                                                                       48028, 23848, 48308,
                                                                       10306, 10414, 26728,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 53188, 0, 3,
                                                                       48308, 24016, 48588,
                                                                       10414, 10522, 26944,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 53548, 0, 3,
                                                                       48588, 24184, 48868,
                                                                       10522, 10630, 27160,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 53908, 0, 3,
                                                                       49148, 24520, 49428,
                                                                       10846, 10954, 27376,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 54268, 0, 3,
                                                                       49428, 24688, 49708,
                                                                       10954, 11062, 27592,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 54628, 0, 3,
                                                                       49708, 24856, 49988,
                                                                       11062, 11170, 27808,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 54988, 0, 3,
                                                                       49988, 25024, 50268,
                                                                       11170, 11278, 28024,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 55348, 0, 3,
                                                                       50268, 25192, 50548,
                                                                       11278, 11386, 28240,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 55708, 0, 3,
                                                                       50548, 25360, 50828,
                                                                       11386, 11494, 28456,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 56068, 0, 3,
                                                                       50828, 25528, 51108,
                                                                       11494, 11602, 28672,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 56428, 0, 3,
                                                                       51388, 25864, 51748,
                                                                       11818, 11953, 28888,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 56878, 0, 3,
                                                                       51748, 26080, 52108,
                                                                       11953, 12088, 29158,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 57328, 0, 3,
                                                                       52108, 26296, 52468,
                                                                       12088, 12223, 29428,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 57778, 0, 3,
                                                                       52468, 26512, 52828,
                                                                       12223, 12358, 29698,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 58228, 0, 3,
                                                                       52828, 26728, 53188,
                                                                       12358, 12493, 29968,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 58678, 0, 3,
                                                                       53188, 26944, 53548,
                                                                       12493, 12628, 30238,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 59128, 0, 3,
                                                                       53908, 27376, 54268,
                                                                       12898, 13033, 30508,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 59578, 0, 3,
                                                                       54268, 27592, 54628,
                                                                       13033, 13168, 30778,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 60028, 0, 3,
                                                                       54628, 27808, 54988,
                                                                       13168, 13303, 31048,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 60478, 0, 3,
                                                                       54988, 28024, 55348,
                                                                       13303, 13438, 31318,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 60928, 0, 3,
                                                                       55348, 28240, 55708,
                                                                       13438, 13573, 31588,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 61378, 0, 3,
                                                                       55708, 28456, 56068,
                                                                       13573, 13708, 31858,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 61828, 0, 3,
                                                                       56428, 28888, 56878,
                                                                       13978, 14143, 32128,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 62378, 0, 3,
                                                                       56878, 29158, 57328,
                                                                       14143, 14308, 32458,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 62928, 0, 3,
                                                                       57328, 29428, 57778,
                                                                       14308, 14473, 32788,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 63478, 0, 3,
                                                                       57778, 29698, 58228,
                                                                       14473, 14638, 33118,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 64028, 0, 3,
                                                                       58228, 29968, 58678,
                                                                       14638, 14803, 33448,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 64578, 0, 3,
                                                                       59128, 30508, 59578,
                                                                       15133, 15298, 33778,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 65128, 0, 3,
                                                                       59578, 30778, 60028,
                                                                       15298, 15463, 34108,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 65678, 0, 3,
                                                                       60028, 31048, 60478,
                                                                       15463, 15628, 34438,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 66228, 0, 3,
                                                                       60478, 31318, 60928,
                                                                       15628, 15793, 34768,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 66778, 0, 3,
                                                                       60928, 31588, 61378,
                                                                       15793, 15958, 35098,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67328, 3, 16288,
                                                                       16294, 35448, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67343, 3, 16294,
                                                                       16300, 35458, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67358, 3, 16300,
                                                                       16306, 35468, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67373, 3, 16306,
                                                                       16312, 35478, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67388, 3, 16312,
                                                                       16318, 35488, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67403, 3, 16318,
                                                                       16324, 35498, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67418, 3, 16324,
                                                                       16330, 35508, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67433, 3, 16330,
                                                                       16336, 35518, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67448, 3, 16336,
                                                                       16342, 35528, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67463, 3, 16342,
                                                                       16348, 35538, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67478, 3, 16348,
                                                                       16354, 35548, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67493, 3, 16354,
                                                                       16360, 35558, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67508, 3, 16372,
                                                                       16378, 35588, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67523, 3, 16378,
                                                                       16384, 35598, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67538, 3, 16384,
                                                                       16390, 35608, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67553, 3, 16390,
                                                                       16396, 35618, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67568, 3, 16396,
                                                                       16402, 35628, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67583, 3, 16402,
                                                                       16408, 35638, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67598, 3, 16408,
                                                                       16414, 35648, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67613, 3, 16414,
                                                                       16420, 35658, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67628, 3, 16420,
                                                                       16426, 35668, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67643, 3, 16426,
                                                                       16432, 35678, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67658, 3, 16432,
                                                                       16438, 35688, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67673, 3, 16438,
                                                                       16444, 35698, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 67688, 0, 3,
                                                                       67328, 35448, 67343,
                                                                       16456, 16474, 35768,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 67733, 0, 3,
                                                                       67343, 35458, 67358,
                                                                       16474, 16492, 35798,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 67778, 0, 3,
                                                                       67358, 35468, 67373,
                                                                       16492, 16510, 35828,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 67823, 0, 3,
                                                                       67373, 35478, 67388,
                                                                       16510, 16528, 35858,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 67868, 0, 3,
                                                                       67388, 35488, 67403,
                                                                       16528, 16546, 35888,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 67913, 0, 3,
                                                                       67403, 35498, 67418,
                                                                       16546, 16564, 35918,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 67958, 0, 3,
                                                                       67418, 35508, 67433,
                                                                       16564, 16582, 35948,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68003, 0, 3,
                                                                       67433, 35518, 67448,
                                                                       16582, 16600, 35978,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68048, 0, 3,
                                                                       67448, 35528, 67463,
                                                                       16600, 16618, 36008,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68093, 0, 3,
                                                                       67463, 35538, 67478,
                                                                       16618, 16636, 36038,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68138, 0, 3,
                                                                       67478, 35548, 67493,
                                                                       16636, 16654, 36068,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68183, 0, 3,
                                                                       67508, 35588, 67523,
                                                                       16690, 16708, 36158,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68228, 0, 3,
                                                                       67523, 35598, 67538,
                                                                       16708, 16726, 36188,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68273, 0, 3,
                                                                       67538, 35608, 67553,
                                                                       16726, 16744, 36218,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68318, 0, 3,
                                                                       67553, 35618, 67568,
                                                                       16744, 16762, 36248,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68363, 0, 3,
                                                                       67568, 35628, 67583,
                                                                       16762, 16780, 36278,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68408, 0, 3,
                                                                       67583, 35638, 67598,
                                                                       16780, 16798, 36308,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68453, 0, 3,
                                                                       67598, 35648, 67613,
                                                                       16798, 16816, 36338,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68498, 0, 3,
                                                                       67613, 35658, 67628,
                                                                       16816, 16834, 36368,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68543, 0, 3,
                                                                       67628, 35668, 67643,
                                                                       16834, 16852, 36398,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68588, 0, 3,
                                                                       67643, 35678, 67658,
                                                                       16852, 16870, 36428,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68633, 0, 3,
                                                                       67658, 35688, 67673,
                                                                       16870, 16888, 36458,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 68678, 0, 3,
                                                                       67688, 35768, 67733,
                                                                       16924, 16960, 36608,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 68768, 0, 3,
                                                                       67733, 35798, 67778,
                                                                       16960, 16996, 36668,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 68858, 0, 3,
                                                                       67778, 35828, 67823,
                                                                       16996, 17032, 36728,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 68948, 0, 3,
                                                                       67823, 35858, 67868,
                                                                       17032, 17068, 36788,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69038, 0, 3,
                                                                       67868, 35888, 67913,
                                                                       17068, 17104, 36848,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69128, 0, 3,
                                                                       67913, 35918, 67958,
                                                                       17104, 17140, 36908,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69218, 0, 3,
                                                                       67958, 35948, 68003,
                                                                       17140, 17176, 36968,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69308, 0, 3,
                                                                       68003, 35978, 68048,
                                                                       17176, 17212, 37028,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69398, 0, 3,
                                                                       68048, 36008, 68093,
                                                                       17212, 17248, 37088,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69488, 0, 3,
                                                                       68093, 36038, 68138,
                                                                       17248, 17284, 37148,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69578, 0, 3,
                                                                       68183, 36158, 68228,
                                                                       17356, 17392, 37328,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69668, 0, 3,
                                                                       68228, 36188, 68273,
                                                                       17392, 17428, 37388,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69758, 0, 3,
                                                                       68273, 36218, 68318,
                                                                       17428, 17464, 37448,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69848, 0, 3,
                                                                       68318, 36248, 68363,
                                                                       17464, 17500, 37508,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69938, 0, 3,
                                                                       68363, 36278, 68408,
                                                                       17500, 17536, 37568,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 70028, 0, 3,
                                                                       68408, 36308, 68453,
                                                                       17536, 17572, 37628,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 70118, 0, 3,
                                                                       68453, 36338, 68498,
                                                                       17572, 17608, 37688,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 70208, 0, 3,
                                                                       68498, 36368, 68543,
                                                                       17608, 17644, 37748,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 70298, 0, 3,
                                                                       68543, 36398, 68588,
                                                                       17644, 17680, 37808,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 70388, 0, 3,
                                                                       68588, 36428, 68633,
                                                                       17680, 17716, 37868,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 70478, 0, 3,
                                                                       68678, 36608, 68768,
                                                                       17788, 17848, 38128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 70628, 0, 3,
                                                                       68768, 36668, 68858,
                                                                       17848, 17908, 38228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 70778, 0, 3,
                                                                       68858, 36728, 68948,
                                                                       17908, 17968, 38328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 70928, 0, 3,
                                                                       68948, 36788, 69038,
                                                                       17968, 18028, 38428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 71078, 0, 3,
                                                                       69038, 36848, 69128,
                                                                       18028, 18088, 38528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 71228, 0, 3,
                                                                       69128, 36908, 69218,
                                                                       18088, 18148, 38628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 71378, 0, 3,
                                                                       69218, 36968, 69308,
                                                                       18148, 18208, 38728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 71528, 0, 3,
                                                                       69308, 37028, 69398,
                                                                       18208, 18268, 38828,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 71678, 0, 3,
                                                                       69398, 37088, 69488,
                                                                       18268, 18328, 38928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 71828, 0, 3,
                                                                       69578, 37328, 69668,
                                                                       18448, 18508, 39228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 71978, 0, 3,
                                                                       69668, 37388, 69758,
                                                                       18508, 18568, 39328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 72128, 0, 3,
                                                                       69758, 37448, 69848,
                                                                       18568, 18628, 39428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 72278, 0, 3,
                                                                       69848, 37508, 69938,
                                                                       18628, 18688, 39528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 72428, 0, 3,
                                                                       69938, 37568, 70028,
                                                                       18688, 18748, 39628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 72578, 0, 3,
                                                                       70028, 37628, 70118,
                                                                       18748, 18808, 39728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 72728, 0, 3,
                                                                       70118, 37688, 70208,
                                                                       18808, 18868, 39828,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 72878, 0, 3,
                                                                       70208, 37748, 70298,
                                                                       18868, 18928, 39928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 73028, 0, 3,
                                                                       70298, 37808, 70388,
                                                                       18928, 18988, 40028,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 73178, 0, 3,
                                                                       70478, 38128, 70628,
                                                                       19108, 19198, 40428,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 73403, 0, 3,
                                                                       70628, 38228, 70778,
                                                                       19198, 19288, 40578,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 73628, 0, 3,
                                                                       70778, 38328, 70928,
                                                                       19288, 19378, 40728,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 73853, 0, 3,
                                                                       70928, 38428, 71078,
                                                                       19378, 19468, 40878,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 74078, 0, 3,
                                                                       71078, 38528, 71228,
                                                                       19468, 19558, 41028,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 74303, 0, 3,
                                                                       71228, 38628, 71378,
                                                                       19558, 19648, 41178,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 74528, 0, 3,
                                                                       71378, 38728, 71528,
                                                                       19648, 19738, 41328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 74753, 0, 3,
                                                                       71528, 38828, 71678,
                                                                       19738, 19828, 41478,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 74978, 0, 3,
                                                                       71828, 39228, 71978,
                                                                       20008, 20098, 41928,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 75203, 0, 3,
                                                                       71978, 39328, 72128,
                                                                       20098, 20188, 42078,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 75428, 0, 3,
                                                                       72128, 39428, 72278,
                                                                       20188, 20278, 42228,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 75653, 0, 3,
                                                                       72278, 39528, 72428,
                                                                       20278, 20368, 42378,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 75878, 0, 3,
                                                                       72428, 39628, 72578,
                                                                       20368, 20458, 42528,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 76103, 0, 3,
                                                                       72578, 39728, 72728,
                                                                       20458, 20548, 42678,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 76328, 0, 3,
                                                                       72728, 39828, 72878,
                                                                       20548, 20638, 42828,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 76553, 0, 3,
                                                                       72878, 39928, 73028,
                                                                       20638, 20728, 42978,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 76778, 0, 3,
                                                                       73178, 40428, 73403,
                                                                       20908, 21034, 43548,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 77093, 0, 3,
                                                                       73403, 40578, 73628,
                                                                       21034, 21160, 43758,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 77408, 0, 3,
                                                                       73628, 40728, 73853,
                                                                       21160, 21286, 43968,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 77723, 0, 3,
                                                                       73853, 40878, 74078,
                                                                       21286, 21412, 44178,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 78038, 0, 3,
                                                                       74078, 41028, 74303,
                                                                       21412, 21538, 44388,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 78353, 0, 3,
                                                                       74303, 41178, 74528,
                                                                       21538, 21664, 44598,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 78668, 0, 3,
                                                                       74528, 41328, 74753,
                                                                       21664, 21790, 44808,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 78983, 0, 3,
                                                                       74978, 41928, 75203,
                                                                       22042, 22168, 45438,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 79298, 0, 3,
                                                                       75203, 42078, 75428,
                                                                       22168, 22294, 45648,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 79613, 0, 3,
                                                                       75428, 42228, 75653,
                                                                       22294, 22420, 45858,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 79928, 0, 3,
                                                                       75653, 42378, 75878,
                                                                       22420, 22546, 46068,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 80243, 0, 3,
                                                                       75878, 42528, 76103,
                                                                       22546, 22672, 46278,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 80558, 0, 3,
                                                                       76103, 42678, 76328,
                                                                       22672, 22798, 46488,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 80873, 0, 3,
                                                                       76328, 42828, 76553,
                                                                       22798, 22924, 46698,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 81188, 0, 3,
                                                                       76778, 43548, 77093,
                                                                       23176, 23344, 47468,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 81608, 0, 3,
                                                                       77093, 43758, 77408,
                                                                       23344, 23512, 47748,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 82028, 0, 3,
                                                                       77408, 43968, 77723,
                                                                       23512, 23680, 48028,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 82448, 0, 3,
                                                                       77723, 44178, 78038,
                                                                       23680, 23848, 48308,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 82868, 0, 3,
                                                                       78038, 44388, 78353,
                                                                       23848, 24016, 48588,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 83288, 0, 3,
                                                                       78353, 44598, 78668,
                                                                       24016, 24184, 48868,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 83708, 0, 3,
                                                                       78983, 45438, 79298,
                                                                       24520, 24688, 49708,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 84128, 0, 3,
                                                                       79298, 45648, 79613,
                                                                       24688, 24856, 49988,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 84548, 0, 3,
                                                                       79613, 45858, 79928,
                                                                       24856, 25024, 50268,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 84968, 0, 3,
                                                                       79928, 46068, 80243,
                                                                       25024, 25192, 50548,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 85388, 0, 3,
                                                                       80243, 46278, 80558,
                                                                       25192, 25360, 50828,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 85808, 0, 3,
                                                                       80558, 46488, 80873,
                                                                       25360, 25528, 51108,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 86228, 0, 3,
                                                                       81188, 47468, 81608,
                                                                       25864, 26080, 52108,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 86768, 0, 3,
                                                                       81608, 47748, 82028,
                                                                       26080, 26296, 52468,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 87308, 0, 3,
                                                                       82028, 48028, 82448,
                                                                       26296, 26512, 52828,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 87848, 0, 3,
                                                                       82448, 48308, 82868,
                                                                       26512, 26728, 53188,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 88388, 0, 3,
                                                                       82868, 48588, 83288,
                                                                       26728, 26944, 53548,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 88928, 0, 3,
                                                                       83708, 49708, 84128,
                                                                       27376, 27592, 54628,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 89468, 0, 3,
                                                                       84128, 49988, 84548,
                                                                       27592, 27808, 54988,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 90008, 0, 3,
                                                                       84548, 50268, 84968,
                                                                       27808, 28024, 55348,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 90548, 0, 3,
                                                                       84968, 50548, 85388,
                                                                       28024, 28240, 55708,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 91088, 0, 3,
                                                                       85388, 50828, 85808,
                                                                       28240, 28456, 56068,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 91628, 0, 3,
                                                                       86228, 52108, 86768,
                                                                       28888, 29158, 57328,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 92303, 0, 3,
                                                                       86768, 52468, 87308,
                                                                       29158, 29428, 57778,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 92978, 0, 3,
                                                                       87308, 52828, 87848,
                                                                       29428, 29698, 58228,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 93653, 0, 3,
                                                                       87848, 53188, 88388,
                                                                       29698, 29968, 58678,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 94328, 0, 3,
                                                                       88928, 54628, 89468,
                                                                       30508, 30778, 60028,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 95003, 0, 3,
                                                                       89468, 54988, 90008,
                                                                       30778, 31048, 60478,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 95678, 0, 3,
                                                                       90008, 55348, 90548,
                                                                       31048, 31318, 60928,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 96353, 0, 3,
                                                                       90548, 55708, 91088,
                                                                       31318, 31588, 61378,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 97028, 0, 3,
                                                                       91628, 57328, 92303,
                                                                       32128, 32458, 62928,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 97853, 0, 3,
                                                                       92303, 57778, 92978,
                                                                       32458, 32788, 63478,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 98678, 0, 3,
                                                                       92978, 58228, 93653,
                                                                       32788, 33118, 64028,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 99503, 0, 3,
                                                                       94328, 60028, 95003,
                                                                       33778, 34108, 65678,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 100328, 0, 3,
                                                                       95003, 60478, 95678,
                                                                       34108, 34438, 66228,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 101153, 0, 3,
                                                                       95678, 60928, 96353,
                                                                       34438, 34768, 66778,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 101978, 3, 35428,
                                                                       35438, 67328, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 101999, 3, 35438,
                                                                       35448, 67343, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102020, 3, 35448,
                                                                       35458, 67358, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102041, 3, 35458,
                                                                       35468, 67373, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102062, 3, 35468,
                                                                       35478, 67388, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102083, 3, 35478,
                                                                       35488, 67403, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102104, 3, 35488,
                                                                       35498, 67418, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102125, 3, 35498,
                                                                       35508, 67433, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102146, 3, 35508,
                                                                       35518, 67448, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102167, 3, 35518,
                                                                       35528, 67463, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102188, 3, 35528,
                                                                       35538, 67478, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102209, 3, 35538,
                                                                       35548, 67493, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102230, 3, 35568,
                                                                       35578, 67508, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102251, 3, 35578,
                                                                       35588, 67523, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102272, 3, 35588,
                                                                       35598, 67538, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102293, 3, 35598,
                                                                       35608, 67553, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102314, 3, 35608,
                                                                       35618, 67568, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102335, 3, 35618,
                                                                       35628, 67583, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102356, 3, 35628,
                                                                       35638, 67598, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102377, 3, 35638,
                                                                       35648, 67613, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102398, 3, 35648,
                                                                       35658, 67628, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102419, 3, 35658,
                                                                       35668, 67643, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102440, 3, 35668,
                                                                       35678, 67658, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102461, 3, 35678,
                                                                       35688, 67673, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 102482, 0, 3,
                                                                       101978, 67328, 101999,
                                                                       35708, 35738, 67688,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 102545, 0, 3,
                                                                       101999, 67343, 102020,
                                                                       35738, 35768, 67733,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 102608, 0, 3,
                                                                       102020, 67358, 102041,
                                                                       35768, 35798, 67778,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 102671, 0, 3,
                                                                       102041, 67373, 102062,
                                                                       35798, 35828, 67823,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 102734, 0, 3,
                                                                       102062, 67388, 102083,
                                                                       35828, 35858, 67868,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 102797, 0, 3,
                                                                       102083, 67403, 102104,
                                                                       35858, 35888, 67913,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 102860, 0, 3,
                                                                       102104, 67418, 102125,
                                                                       35888, 35918, 67958,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 102923, 0, 3,
                                                                       102125, 67433, 102146,
                                                                       35918, 35948, 68003,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 102986, 0, 3,
                                                                       102146, 67448, 102167,
                                                                       35948, 35978, 68048,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 103049, 0, 3,
                                                                       102167, 67463, 102188,
                                                                       35978, 36008, 68093,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 103112, 0, 3,
                                                                       102188, 67478, 102209,
                                                                       36008, 36038, 68138,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 103175, 0, 3,
                                                                       102230, 67508, 102251,
                                                                       36098, 36128, 68183,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 103238, 0, 3,
                                                                       102251, 67523, 102272,
                                                                       36128, 36158, 68228,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 103301, 0, 3,
                                                                       102272, 67538, 102293,
                                                                       36158, 36188, 68273,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 103364, 0, 3,
                                                                       102293, 67553, 102314,
                                                                       36188, 36218, 68318,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 103427, 0, 3,
                                                                       102314, 67568, 102335,
                                                                       36218, 36248, 68363,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 103490, 0, 3,
                                                                       102335, 67583, 102356,
                                                                       36248, 36278, 68408,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 103553, 0, 3,
                                                                       102356, 67598, 102377,
                                                                       36278, 36308, 68453,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 103616, 0, 3,
                                                                       102377, 67613, 102398,
                                                                       36308, 36338, 68498,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 103679, 0, 3,
                                                                       102398, 67628, 102419,
                                                                       36338, 36368, 68543,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 103742, 0, 3,
                                                                       102419, 67643, 102440,
                                                                       36368, 36398, 68588,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 103805, 0, 3,
                                                                       102440, 67658, 102461,
                                                                       36398, 36428, 68633,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 103868, 0, 3,
                                                                       102482, 67688, 102545,
                                                                       36488, 36548, 68678,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 103994, 0, 3,
                                                                       102545, 67733, 102608,
                                                                       36548, 36608, 68768,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 104120, 0, 3,
                                                                       102608, 67778, 102671,
                                                                       36608, 36668, 68858,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 104246, 0, 3,
                                                                       102671, 67823, 102734,
                                                                       36668, 36728, 68948,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 104372, 0, 3,
                                                                       102734, 67868, 102797,
                                                                       36728, 36788, 69038,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 104498, 0, 3,
                                                                       102797, 67913, 102860,
                                                                       36788, 36848, 69128,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 104624, 0, 3,
                                                                       102860, 67958, 102923,
                                                                       36848, 36908, 69218,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 104750, 0, 3,
                                                                       102923, 68003, 102986,
                                                                       36908, 36968, 69308,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 104876, 0, 3,
                                                                       102986, 68048, 103049,
                                                                       36968, 37028, 69398,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 105002, 0, 3,
                                                                       103049, 68093, 103112,
                                                                       37028, 37088, 69488,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 105128, 0, 3,
                                                                       103175, 68183, 103238,
                                                                       37208, 37268, 69578,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 105254, 0, 3,
                                                                       103238, 68228, 103301,
                                                                       37268, 37328, 69668,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 105380, 0, 3,
                                                                       103301, 68273, 103364,
                                                                       37328, 37388, 69758,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 105506, 0, 3,
                                                                       103364, 68318, 103427,
                                                                       37388, 37448, 69848,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 105632, 0, 3,
                                                                       103427, 68363, 103490,
                                                                       37448, 37508, 69938,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 105758, 0, 3,
                                                                       103490, 68408, 103553,
                                                                       37508, 37568, 70028,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 105884, 0, 3,
                                                                       103553, 68453, 103616,
                                                                       37568, 37628, 70118,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 106010, 0, 3,
                                                                       103616, 68498, 103679,
                                                                       37628, 37688, 70208,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 106136, 0, 3,
                                                                       103679, 68543, 103742,
                                                                       37688, 37748, 70298,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 106262, 0, 3,
                                                                       103742, 68588, 103805,
                                                                       37748, 37808, 70388,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 106388, 0, 3,
                                                                       103868, 68678, 103994,
                                                                       37928, 38028, 70478,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 106598, 0, 3,
                                                                       103994, 68768, 104120,
                                                                       38028, 38128, 70628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 106808, 0, 3,
                                                                       104120, 68858, 104246,
                                                                       38128, 38228, 70778,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 107018, 0, 3,
                                                                       104246, 68948, 104372,
                                                                       38228, 38328, 70928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 107228, 0, 3,
                                                                       104372, 69038, 104498,
                                                                       38328, 38428, 71078,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 107438, 0, 3,
                                                                       104498, 69128, 104624,
                                                                       38428, 38528, 71228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 107648, 0, 3,
                                                                       104624, 69218, 104750,
                                                                       38528, 38628, 71378,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 107858, 0, 3,
                                                                       104750, 69308, 104876,
                                                                       38628, 38728, 71528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 108068, 0, 3,
                                                                       104876, 69398, 105002,
                                                                       38728, 38828, 71678,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 108278, 0, 3,
                                                                       105128, 69578, 105254,
                                                                       39028, 39128, 71828,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 108488, 0, 3,
                                                                       105254, 69668, 105380,
                                                                       39128, 39228, 71978,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 108698, 0, 3,
                                                                       105380, 69758, 105506,
                                                                       39228, 39328, 72128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 108908, 0, 3,
                                                                       105506, 69848, 105632,
                                                                       39328, 39428, 72278,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 109118, 0, 3,
                                                                       105632, 69938, 105758,
                                                                       39428, 39528, 72428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 109328, 0, 3,
                                                                       105758, 70028, 105884,
                                                                       39528, 39628, 72578,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 109538, 0, 3,
                                                                       105884, 70118, 106010,
                                                                       39628, 39728, 72728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 109748, 0, 3,
                                                                       106010, 70208, 106136,
                                                                       39728, 39828, 72878,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 109958, 0, 3,
                                                                       106136, 70298, 106262,
                                                                       39828, 39928, 73028,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 110168, 0, 3,
                                                                       106388, 70478, 106598,
                                                                       40128, 40278, 73178,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 110483, 0, 3,
                                                                       106598, 70628, 106808,
                                                                       40278, 40428, 73403,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 110798, 0, 3,
                                                                       106808, 70778, 107018,
                                                                       40428, 40578, 73628,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 111113, 0, 3,
                                                                       107018, 70928, 107228,
                                                                       40578, 40728, 73853,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 111428, 0, 3,
                                                                       107228, 71078, 107438,
                                                                       40728, 40878, 74078,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 111743, 0, 3,
                                                                       107438, 71228, 107648,
                                                                       40878, 41028, 74303,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 112058, 0, 3,
                                                                       107648, 71378, 107858,
                                                                       41028, 41178, 74528,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 112373, 0, 3,
                                                                       107858, 71528, 108068,
                                                                       41178, 41328, 74753,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 112688, 0, 3,
                                                                       108278, 71828, 108488,
                                                                       41628, 41778, 74978,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 113003, 0, 3,
                                                                       108488, 71978, 108698,
                                                                       41778, 41928, 75203,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 113318, 0, 3,
                                                                       108698, 72128, 108908,
                                                                       41928, 42078, 75428,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 113633, 0, 3,
                                                                       108908, 72278, 109118,
                                                                       42078, 42228, 75653,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 113948, 0, 3,
                                                                       109118, 72428, 109328,
                                                                       42228, 42378, 75878,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 114263, 0, 3,
                                                                       109328, 72578, 109538,
                                                                       42378, 42528, 76103,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 114578, 0, 3,
                                                                       109538, 72728, 109748,
                                                                       42528, 42678, 76328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 114893, 0, 3,
                                                                       109748, 72878, 109958,
                                                                       42678, 42828, 76553,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 115208, 0, 3,
                                                                       110168, 73178, 110483,
                                                                       43128, 43338, 76778,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 115649, 0, 3,
                                                                       110483, 73403, 110798,
                                                                       43338, 43548, 77093,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 116090, 0, 3,
                                                                       110798, 73628, 111113,
                                                                       43548, 43758, 77408,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 116531, 0, 3,
                                                                       111113, 73853, 111428,
                                                                       43758, 43968, 77723,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 116972, 0, 3,
                                                                       111428, 74078, 111743,
                                                                       43968, 44178, 78038,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 117413, 0, 3,
                                                                       111743, 74303, 112058,
                                                                       44178, 44388, 78353,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 117854, 0, 3,
                                                                       112058, 74528, 112373,
                                                                       44388, 44598, 78668,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 118295, 0, 3,
                                                                       112688, 74978, 113003,
                                                                       45018, 45228, 78983,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 118736, 0, 3,
                                                                       113003, 75203, 113318,
                                                                       45228, 45438, 79298,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 119177, 0, 3,
                                                                       113318, 75428, 113633,
                                                                       45438, 45648, 79613,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 119618, 0, 3,
                                                                       113633, 75653, 113948,
                                                                       45648, 45858, 79928,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 120059, 0, 3,
                                                                       113948, 75878, 114263,
                                                                       45858, 46068, 80243,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 120500, 0, 3,
                                                                       114263, 76103, 114578,
                                                                       46068, 46278, 80558,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 120941, 0, 3,
                                                                       114578, 76328, 114893,
                                                                       46278, 46488, 80873,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 121382, 0, 3,
                                                                       115208, 76778, 115649,
                                                                       46908, 47188, 81188,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 121970, 0, 3,
                                                                       115649, 77093, 116090,
                                                                       47188, 47468, 81608,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 122558, 0, 3,
                                                                       116090, 77408, 116531,
                                                                       47468, 47748, 82028,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 123146, 0, 3,
                                                                       116531, 77723, 116972,
                                                                       47748, 48028, 82448,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 123734, 0, 3,
                                                                       116972, 78038, 117413,
                                                                       48028, 48308, 82868,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 124322, 0, 3,
                                                                       117413, 78353, 117854,
                                                                       48308, 48588, 83288,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 124910, 0, 3,
                                                                       118295, 78983, 118736,
                                                                       49148, 49428, 83708,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 125498, 0, 3,
                                                                       118736, 79298, 119177,
                                                                       49428, 49708, 84128,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 126086, 0, 3,
                                                                       119177, 79613, 119618,
                                                                       49708, 49988, 84548,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 126674, 0, 3,
                                                                       119618, 79928, 120059,
                                                                       49988, 50268, 84968,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 127262, 0, 3,
                                                                       120059, 80243, 120500,
                                                                       50268, 50548, 85388,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 127850, 0, 3,
                                                                       120500, 80558, 120941,
                                                                       50548, 50828, 85808,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 128438, 0, 3,
                                                                       121382, 81188, 121970,
                                                                       51388, 51748, 86228,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 129194, 0, 3,
                                                                       121970, 81608, 122558,
                                                                       51748, 52108, 86768,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 129950, 0, 3,
                                                                       122558, 82028, 123146,
                                                                       52108, 52468, 87308,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 130706, 0, 3,
                                                                       123146, 82448, 123734,
                                                                       52468, 52828, 87848,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 131462, 0, 3,
                                                                       123734, 82868, 124322,
                                                                       52828, 53188, 88388,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 132218, 0, 3,
                                                                       124910, 83708, 125498,
                                                                       53908, 54268, 88928,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 132974, 0, 3,
                                                                       125498, 84128, 126086,
                                                                       54268, 54628, 89468,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 133730, 0, 3,
                                                                       126086, 84548, 126674,
                                                                       54628, 54988, 90008,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 134486, 0, 3,
                                                                       126674, 84968, 127262,
                                                                       54988, 55348, 90548,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 135242, 0, 3,
                                                                       127262, 85388, 127850,
                                                                       55348, 55708, 91088,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 135998, 0, 3,
                                                                       128438, 86228, 129194,
                                                                       56428, 56878, 91628,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 136943, 0, 3,
                                                                       129194, 86768, 129950,
                                                                       56878, 57328, 92303,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 137888, 0, 3,
                                                                       129950, 87308, 130706,
                                                                       57328, 57778, 92978,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 138833, 0, 3,
                                                                       130706, 87848, 131462,
                                                                       57778, 58228, 93653,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 139778, 0, 3,
                                                                       132218, 88928, 132974,
                                                                       59128, 59578, 94328,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 140723, 0, 3,
                                                                       132974, 89468, 133730,
                                                                       59578, 60028, 95003,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 141668, 0, 3,
                                                                       133730, 90008, 134486,
                                                                       60028, 60478, 95678,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 142613, 0, 3,
                                                                       134486, 90548, 135242,
                                                                       60478, 60928, 96353,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 143558, 0, 3,
                                                                       135998, 91628, 136943,
                                                                       61828, 62378, 97028,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 144713, 0, 3,
                                                                       136943, 92303, 137888,
                                                                       62378, 62928, 97853,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 145868, 0, 3,
                                                                       137888, 92978, 138833,
                                                                       62928, 63478, 98678,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 147023, 0, 3,
                                                                       139778, 94328, 140723,
                                                                       64578, 65128, 99503,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 148178, 0, 3,
                                                                       140723, 95003, 141668,
                                                                       65128, 65678, 100328,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 149333, 0, 3,
                                                                       141668, 95678, 142613,
                                                                       65678, 66228, 101153,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150488, 3, 67328,
                                                                       67343, 102020, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150516, 3, 67343,
                                                                       67358, 102041, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150544, 3, 67358,
                                                                       67373, 102062, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150572, 3, 67373,
                                                                       67388, 102083, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150600, 3, 67388,
                                                                       67403, 102104, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150628, 3, 67403,
                                                                       67418, 102125, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150656, 3, 67418,
                                                                       67433, 102146, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150684, 3, 67433,
                                                                       67448, 102167, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150712, 3, 67448,
                                                                       67463, 102188, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150740, 3, 67463,
                                                                       67478, 102209, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150768, 3, 67508,
                                                                       67523, 102272, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150796, 3, 67523,
                                                                       67538, 102293, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150824, 3, 67538,
                                                                       67553, 102314, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150852, 3, 67553,
                                                                       67568, 102335, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150880, 3, 67568,
                                                                       67583, 102356, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150908, 3, 67583,
                                                                       67598, 102377, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150936, 3, 67598,
                                                                       67613, 102398, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150964, 3, 67613,
                                                                       67628, 102419, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150992, 3, 67628,
                                                                       67643, 102440, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 151020, 3, 67643,
                                                                       67658, 102461, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 151048, 0, 3,
                                                                       150488, 102020, 150516,
                                                                       67688, 67733, 102608,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 151132, 0, 3,
                                                                       150516, 102041, 150544,
                                                                       67733, 67778, 102671,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 151216, 0, 3,
                                                                       150544, 102062, 150572,
                                                                       67778, 67823, 102734,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 151300, 0, 3,
                                                                       150572, 102083, 150600,
                                                                       67823, 67868, 102797,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 151384, 0, 3,
                                                                       150600, 102104, 150628,
                                                                       67868, 67913, 102860,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 151468, 0, 3,
                                                                       150628, 102125, 150656,
                                                                       67913, 67958, 102923,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 151552, 0, 3,
                                                                       150656, 102146, 150684,
                                                                       67958, 68003, 102986,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 151636, 0, 3,
                                                                       150684, 102167, 150712,
                                                                       68003, 68048, 103049,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 151720, 0, 3,
                                                                       150712, 102188, 150740,
                                                                       68048, 68093, 103112,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 151804, 0, 3,
                                                                       150768, 102272, 150796,
                                                                       68183, 68228, 103301,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 151888, 0, 3,
                                                                       150796, 102293, 150824,
                                                                       68228, 68273, 103364,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 151972, 0, 3,
                                                                       150824, 102314, 150852,
                                                                       68273, 68318, 103427,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 152056, 0, 3,
                                                                       150852, 102335, 150880,
                                                                       68318, 68363, 103490,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 152140, 0, 3,
                                                                       150880, 102356, 150908,
                                                                       68363, 68408, 103553,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 152224, 0, 3,
                                                                       150908, 102377, 150936,
                                                                       68408, 68453, 103616,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 152308, 0, 3,
                                                                       150936, 102398, 150964,
                                                                       68453, 68498, 103679,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 152392, 0, 3,
                                                                       150964, 102419, 150992,
                                                                       68498, 68543, 103742,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 152476, 0, 3,
                                                                       150992, 102440, 151020,
                                                                       68543, 68588, 103805,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 152560, 0, 3,
                                                                       151048, 102608, 151132,
                                                                       68678, 68768, 104120,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 152728, 0, 3,
                                                                       151132, 102671, 151216,
                                                                       68768, 68858, 104246,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 152896, 0, 3,
                                                                       151216, 102734, 151300,
                                                                       68858, 68948, 104372,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 153064, 0, 3,
                                                                       151300, 102797, 151384,
                                                                       68948, 69038, 104498,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 153232, 0, 3,
                                                                       151384, 102860, 151468,
                                                                       69038, 69128, 104624,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 153400, 0, 3,
                                                                       151468, 102923, 151552,
                                                                       69128, 69218, 104750,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 153568, 0, 3,
                                                                       151552, 102986, 151636,
                                                                       69218, 69308, 104876,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 153736, 0, 3,
                                                                       151636, 103049, 151720,
                                                                       69308, 69398, 105002,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 153904, 0, 3,
                                                                       151804, 103301, 151888,
                                                                       69578, 69668, 105380,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 154072, 0, 3,
                                                                       151888, 103364, 151972,
                                                                       69668, 69758, 105506,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 154240, 0, 3,
                                                                       151972, 103427, 152056,
                                                                       69758, 69848, 105632,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 154408, 0, 3,
                                                                       152056, 103490, 152140,
                                                                       69848, 69938, 105758,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 154576, 0, 3,
                                                                       152140, 103553, 152224,
                                                                       69938, 70028, 105884,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 154744, 0, 3,
                                                                       152224, 103616, 152308,
                                                                       70028, 70118, 106010,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 154912, 0, 3,
                                                                       152308, 103679, 152392,
                                                                       70118, 70208, 106136,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 155080, 0, 3,
                                                                       152392, 103742, 152476,
                                                                       70208, 70298, 106262,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 155248, 0, 3,
                                                                       152560, 104120, 152728,
                                                                       70478, 70628, 106808,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 155528, 0, 3,
                                                                       152728, 104246, 152896,
                                                                       70628, 70778, 107018,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 155808, 0, 3,
                                                                       152896, 104372, 153064,
                                                                       70778, 70928, 107228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 156088, 0, 3,
                                                                       153064, 104498, 153232,
                                                                       70928, 71078, 107438,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 156368, 0, 3,
                                                                       153232, 104624, 153400,
                                                                       71078, 71228, 107648,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 156648, 0, 3,
                                                                       153400, 104750, 153568,
                                                                       71228, 71378, 107858,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 156928, 0, 3,
                                                                       153568, 104876, 153736,
                                                                       71378, 71528, 108068,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 157208, 0, 3,
                                                                       153904, 105380, 154072,
                                                                       71828, 71978, 108698,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 157488, 0, 3,
                                                                       154072, 105506, 154240,
                                                                       71978, 72128, 108908,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 157768, 0, 3,
                                                                       154240, 105632, 154408,
                                                                       72128, 72278, 109118,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 158048, 0, 3,
                                                                       154408, 105758, 154576,
                                                                       72278, 72428, 109328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 158328, 0, 3,
                                                                       154576, 105884, 154744,
                                                                       72428, 72578, 109538,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 158608, 0, 3,
                                                                       154744, 106010, 154912,
                                                                       72578, 72728, 109748,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 158888, 0, 3,
                                                                       154912, 106136, 155080,
                                                                       72728, 72878, 109958,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 159168, 0, 3,
                                                                       155248, 106808, 155528,
                                                                       73178, 73403, 110798,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 159588, 0, 3,
                                                                       155528, 107018, 155808,
                                                                       73403, 73628, 111113,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 160008, 0, 3,
                                                                       155808, 107228, 156088,
                                                                       73628, 73853, 111428,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 160428, 0, 3,
                                                                       156088, 107438, 156368,
                                                                       73853, 74078, 111743,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 160848, 0, 3,
                                                                       156368, 107648, 156648,
                                                                       74078, 74303, 112058,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 161268, 0, 3,
                                                                       156648, 107858, 156928,
                                                                       74303, 74528, 112373,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 161688, 0, 3,
                                                                       157208, 108698, 157488,
                                                                       74978, 75203, 113318,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 162108, 0, 3,
                                                                       157488, 108908, 157768,
                                                                       75203, 75428, 113633,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 162528, 0, 3,
                                                                       157768, 109118, 158048,
                                                                       75428, 75653, 113948,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 162948, 0, 3,
                                                                       158048, 109328, 158328,
                                                                       75653, 75878, 114263,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 163368, 0, 3,
                                                                       158328, 109538, 158608,
                                                                       75878, 76103, 114578,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 163788, 0, 3,
                                                                       158608, 109748, 158888,
                                                                       76103, 76328, 114893,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 164208, 0, 3,
                                                                       159168, 110798, 159588,
                                                                       76778, 77093, 116090,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 164796, 0, 3,
                                                                       159588, 111113, 160008,
                                                                       77093, 77408, 116531,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 165384, 0, 3,
                                                                       160008, 111428, 160428,
                                                                       77408, 77723, 116972,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 165972, 0, 3,
                                                                       160428, 111743, 160848,
                                                                       77723, 78038, 117413,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 166560, 0, 3,
                                                                       160848, 112058, 161268,
                                                                       78038, 78353, 117854,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 167148, 0, 3,
                                                                       161688, 113318, 162108,
                                                                       78983, 79298, 119177,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 167736, 0, 3,
                                                                       162108, 113633, 162528,
                                                                       79298, 79613, 119618,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 168324, 0, 3,
                                                                       162528, 113948, 162948,
                                                                       79613, 79928, 120059,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 168912, 0, 3,
                                                                       162948, 114263, 163368,
                                                                       79928, 80243, 120500,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 169500, 0, 3,
                                                                       163368, 114578, 163788,
                                                                       80243, 80558, 120941,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 170088, 0, 3,
                                                                       164208, 116090, 164796,
                                                                       81188, 81608, 122558,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 170872, 0, 3,
                                                                       164796, 116531, 165384,
                                                                       81608, 82028, 123146,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 171656, 0, 3,
                                                                       165384, 116972, 165972,
                                                                       82028, 82448, 123734,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 172440, 0, 3,
                                                                       165972, 117413, 166560,
                                                                       82448, 82868, 124322,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 173224, 0, 3,
                                                                       167148, 119177, 167736,
                                                                       83708, 84128, 126086,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 174008, 0, 3,
                                                                       167736, 119618, 168324,
                                                                       84128, 84548, 126674,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 174792, 0, 3,
                                                                       168324, 120059, 168912,
                                                                       84548, 84968, 127262,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 175576, 0, 3,
                                                                       168912, 120500, 169500,
                                                                       84968, 85388, 127850,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 176360, 0, 3,
                                                                       170088, 122558, 170872,
                                                                       86228, 86768, 129950,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 177368, 0, 3,
                                                                       170872, 123146, 171656,
                                                                       86768, 87308, 130706,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 178376, 0, 3,
                                                                       171656, 123734, 172440,
                                                                       87308, 87848, 131462,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 179384, 0, 3,
                                                                       173224, 126086, 174008,
                                                                       88928, 89468, 133730,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 180392, 0, 3,
                                                                       174008, 126674, 174792,
                                                                       89468, 90008, 134486,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 181400, 0, 3,
                                                                       174792, 127262, 175576,
                                                                       90008, 90548, 135242,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 182408, 0, 3,
                                                                       176360, 129950, 177368,
                                                                       91628, 92303, 137888,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 183668, 0, 3,
                                                                       177368, 130706, 178376,
                                                                       92303, 92978, 138833,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 184928, 0, 3,
                                                                       179384, 133730, 180392,
                                                                       94328, 95003, 141668,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 186188, 0, 3,
                                                                       180392, 134486, 181400,
                                                                       95003, 95678, 142613,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 187448, 0, 3,
                                                                       182408, 137888, 183668,
                                                                       97028, 97853, 145868,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 188988, 0, 3,
                                                                       184928, 141668, 186188,
                                                                       99503, 100328, 149333,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190528, 3, 101978,
                                                                       101999, 150488, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190564, 3, 101999,
                                                                       102020, 150516, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190600, 3, 102020,
                                                                       102041, 150544, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190636, 3, 102041,
                                                                       102062, 150572, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190672, 3, 102062,
                                                                       102083, 150600, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190708, 3, 102083,
                                                                       102104, 150628, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190744, 3, 102104,
                                                                       102125, 150656, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190780, 3, 102125,
                                                                       102146, 150684, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190816, 3, 102146,
                                                                       102167, 150712, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190852, 3, 102167,
                                                                       102188, 150740, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190888, 3, 102230,
                                                                       102251, 150768, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190924, 3, 102251,
                                                                       102272, 150796, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190960, 3, 102272,
                                                                       102293, 150824, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190996, 3, 102293,
                                                                       102314, 150852, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 191032, 3, 102314,
                                                                       102335, 150880, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 191068, 3, 102335,
                                                                       102356, 150908, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 191104, 3, 102356,
                                                                       102377, 150936, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 191140, 3, 102377,
                                                                       102398, 150964, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 191176, 3, 102398,
                                                                       102419, 150992, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 191212, 3, 102419,
                                                                       102440, 151020, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 191248, 0, 3,
                                                                       190528, 150488, 190564,
                                                                       102482, 102545, 151048,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 191356, 0, 3,
                                                                       190564, 150516, 190600,
                                                                       102545, 102608, 151132,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 191464, 0, 3,
                                                                       190600, 150544, 190636,
                                                                       102608, 102671, 151216,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 191572, 0, 3,
                                                                       190636, 150572, 190672,
                                                                       102671, 102734, 151300,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 191680, 0, 3,
                                                                       190672, 150600, 190708,
                                                                       102734, 102797, 151384,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 191788, 0, 3,
                                                                       190708, 150628, 190744,
                                                                       102797, 102860, 151468,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 191896, 0, 3,
                                                                       190744, 150656, 190780,
                                                                       102860, 102923, 151552,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 192004, 0, 3,
                                                                       190780, 150684, 190816,
                                                                       102923, 102986, 151636,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 192112, 0, 3,
                                                                       190816, 150712, 190852,
                                                                       102986, 103049, 151720,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 192220, 0, 3,
                                                                       190888, 150768, 190924,
                                                                       103175, 103238, 151804,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 192328, 0, 3,
                                                                       190924, 150796, 190960,
                                                                       103238, 103301, 151888,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 192436, 0, 3,
                                                                       190960, 150824, 190996,
                                                                       103301, 103364, 151972,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 192544, 0, 3,
                                                                       190996, 150852, 191032,
                                                                       103364, 103427, 152056,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 192652, 0, 3,
                                                                       191032, 150880, 191068,
                                                                       103427, 103490, 152140,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 192760, 0, 3,
                                                                       191068, 150908, 191104,
                                                                       103490, 103553, 152224,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 192868, 0, 3,
                                                                       191104, 150936, 191140,
                                                                       103553, 103616, 152308,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 192976, 0, 3,
                                                                       191140, 150964, 191176,
                                                                       103616, 103679, 152392,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 193084, 0, 3,
                                                                       191176, 150992, 191212,
                                                                       103679, 103742, 152476,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 193192, 0, 3,
                                                                       191248, 151048, 191356,
                                                                       103868, 103994, 152560,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 193408, 0, 3,
                                                                       191356, 151132, 191464,
                                                                       103994, 104120, 152728,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 193624, 0, 3,
                                                                       191464, 151216, 191572,
                                                                       104120, 104246, 152896,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 193840, 0, 3,
                                                                       191572, 151300, 191680,
                                                                       104246, 104372, 153064,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 194056, 0, 3,
                                                                       191680, 151384, 191788,
                                                                       104372, 104498, 153232,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 194272, 0, 3,
                                                                       191788, 151468, 191896,
                                                                       104498, 104624, 153400,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 194488, 0, 3,
                                                                       191896, 151552, 192004,
                                                                       104624, 104750, 153568,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 194704, 0, 3,
                                                                       192004, 151636, 192112,
                                                                       104750, 104876, 153736,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 194920, 0, 3,
                                                                       192220, 151804, 192328,
                                                                       105128, 105254, 153904,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 195136, 0, 3,
                                                                       192328, 151888, 192436,
                                                                       105254, 105380, 154072,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 195352, 0, 3,
                                                                       192436, 151972, 192544,
                                                                       105380, 105506, 154240,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 195568, 0, 3,
                                                                       192544, 152056, 192652,
                                                                       105506, 105632, 154408,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 195784, 0, 3,
                                                                       192652, 152140, 192760,
                                                                       105632, 105758, 154576,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 196000, 0, 3,
                                                                       192760, 152224, 192868,
                                                                       105758, 105884, 154744,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 196216, 0, 3,
                                                                       192868, 152308, 192976,
                                                                       105884, 106010, 154912,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 196432, 0, 3,
                                                                       192976, 152392, 193084,
                                                                       106010, 106136, 155080,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 196648, 0, 3,
                                                                       193192, 152560, 193408,
                                                                       106388, 106598, 155248,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 197008, 0, 3,
                                                                       193408, 152728, 193624,
                                                                       106598, 106808, 155528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 197368, 0, 3,
                                                                       193624, 152896, 193840,
                                                                       106808, 107018, 155808,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 197728, 0, 3,
                                                                       193840, 153064, 194056,
                                                                       107018, 107228, 156088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 198088, 0, 3,
                                                                       194056, 153232, 194272,
                                                                       107228, 107438, 156368,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 198448, 0, 3,
                                                                       194272, 153400, 194488,
                                                                       107438, 107648, 156648,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 198808, 0, 3,
                                                                       194488, 153568, 194704,
                                                                       107648, 107858, 156928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 199168, 0, 3,
                                                                       194920, 153904, 195136,
                                                                       108278, 108488, 157208,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 199528, 0, 3,
                                                                       195136, 154072, 195352,
                                                                       108488, 108698, 157488,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 199888, 0, 3,
                                                                       195352, 154240, 195568,
                                                                       108698, 108908, 157768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 200248, 0, 3,
                                                                       195568, 154408, 195784,
                                                                       108908, 109118, 158048,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 200608, 0, 3,
                                                                       195784, 154576, 196000,
                                                                       109118, 109328, 158328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 200968, 0, 3,
                                                                       196000, 154744, 196216,
                                                                       109328, 109538, 158608,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 201328, 0, 3,
                                                                       196216, 154912, 196432,
                                                                       109538, 109748, 158888,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 201688, 0, 3,
                                                                       196648, 155248, 197008,
                                                                       110168, 110483, 159168,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 202228, 0, 3,
                                                                       197008, 155528, 197368,
                                                                       110483, 110798, 159588,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 202768, 0, 3,
                                                                       197368, 155808, 197728,
                                                                       110798, 111113, 160008,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 203308, 0, 3,
                                                                       197728, 156088, 198088,
                                                                       111113, 111428, 160428,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 203848, 0, 3,
                                                                       198088, 156368, 198448,
                                                                       111428, 111743, 160848,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 204388, 0, 3,
                                                                       198448, 156648, 198808,
                                                                       111743, 112058, 161268,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 204928, 0, 3,
                                                                       199168, 157208, 199528,
                                                                       112688, 113003, 161688,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 205468, 0, 3,
                                                                       199528, 157488, 199888,
                                                                       113003, 113318, 162108,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 206008, 0, 3,
                                                                       199888, 157768, 200248,
                                                                       113318, 113633, 162528,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 206548, 0, 3,
                                                                       200248, 158048, 200608,
                                                                       113633, 113948, 162948,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 207088, 0, 3,
                                                                       200608, 158328, 200968,
                                                                       113948, 114263, 163368,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 207628, 0, 3,
                                                                       200968, 158608, 201328,
                                                                       114263, 114578, 163788,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 208168, 0, 3,
                                                                       201688, 159168, 202228,
                                                                       115208, 115649, 164208,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 208924, 0, 3,
                                                                       202228, 159588, 202768,
                                                                       115649, 116090, 164796,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 209680, 0, 3,
                                                                       202768, 160008, 203308,
                                                                       116090, 116531, 165384,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 210436, 0, 3,
                                                                       203308, 160428, 203848,
                                                                       116531, 116972, 165972,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 211192, 0, 3,
                                                                       203848, 160848, 204388,
                                                                       116972, 117413, 166560,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 211948, 0, 3,
                                                                       204928, 161688, 205468,
                                                                       118295, 118736, 167148,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 212704, 0, 3,
                                                                       205468, 162108, 206008,
                                                                       118736, 119177, 167736,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 213460, 0, 3,
                                                                       206008, 162528, 206548,
                                                                       119177, 119618, 168324,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 214216, 0, 3,
                                                                       206548, 162948, 207088,
                                                                       119618, 120059, 168912,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 214972, 0, 3,
                                                                       207088, 163368, 207628,
                                                                       120059, 120500, 169500,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 215728, 0, 3,
                                                                       208168, 164208, 208924,
                                                                       121382, 121970, 170088,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 216736, 0, 3,
                                                                       208924, 164796, 209680,
                                                                       121970, 122558, 170872,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 217744, 0, 3,
                                                                       209680, 165384, 210436,
                                                                       122558, 123146, 171656,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 218752, 0, 3,
                                                                       210436, 165972, 211192,
                                                                       123146, 123734, 172440,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 219760, 0, 3,
                                                                       211948, 167148, 212704,
                                                                       124910, 125498, 173224,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 220768, 0, 3,
                                                                       212704, 167736, 213460,
                                                                       125498, 126086, 174008,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 221776, 0, 3,
                                                                       213460, 168324, 214216,
                                                                       126086, 126674, 174792,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 222784, 0, 3,
                                                                       214216, 168912, 214972,
                                                                       126674, 127262, 175576,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 223792, 0, 3,
                                                                       215728, 170088, 216736,
                                                                       128438, 129194, 176360,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 225088, 0, 3,
                                                                       216736, 170872, 217744,
                                                                       129194, 129950, 177368,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 226384, 0, 3,
                                                                       217744, 171656, 218752,
                                                                       129950, 130706, 178376,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 227680, 0, 3,
                                                                       219760, 173224, 220768,
                                                                       132218, 132974, 179384,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 228976, 0, 3,
                                                                       220768, 174008, 221776,
                                                                       132974, 133730, 180392,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 230272, 0, 3,
                                                                       221776, 174792, 222784,
                                                                       133730, 134486, 181400,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 231568, 0, 3,
                                                                       223792, 176360, 225088,
                                                                       135998, 136943, 182408,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 233188, 0, 3,
                                                                       225088, 177368, 226384,
                                                                       136943, 137888, 183668,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 234808, 0, 3,
                                                                       227680, 179384, 228976,
                                                                       139778, 140723, 184928,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 236428, 0, 3,
                                                                       228976, 180392, 230272,
                                                                       140723, 141668, 186188,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 238048, 0, 3,
                                                                       231568, 182408, 233188,
                                                                       143558, 144713, 187448,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 240028, 0, 3,
                                                                       234808, 184928, 236428,
                                                                       147023, 148178, 188988,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 242008, 208168, 756, ncols);

                    simdfunc::contract_primitives(buffer, 243079, 211948, 756, ncols);

                    simdfunc::contract_primitives(buffer, 244150, 215728, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 245578, 219760, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 247006, 223792, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 248842, 227680, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 250678, 231568, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 252973, 234808, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 255268, 238048, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 258073, 240028, 1980, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 242764, 242008, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 243835, 243079, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 245158, 244150, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 246586, 245578, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 248302, 247006, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 250138, 248842, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 252298, 250678, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 254593, 252973, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 257248, 255268, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 260053, 258073, 55, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 260878, 242764, 245158, 15,
                                             nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 261823, 243835, 246586, 15,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 262768, 245158, 248302, 15,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 264028, 246586, 250138, 15,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 265288, 248302, 252298, 15,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 266908, 250138, 254593, 15,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 268528, 252298, 257248, 15,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 270553, 254593, 260053, 15,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 272578, 260878, 262768, 15,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 274468, 261823, 264028, 15,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 276358, 262768, 265288, 15,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 278878, 264028, 266908, 15,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 281398, 265288, 268528, 15,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 284638, 266908, 270553, 15,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 287878, 272578, 276358, 15,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 291028, 274468, 278878, 15,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 294178, 276358, 281398, 15,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 298378, 278878, 284638, 15,
                                             nmax);

        simdtrf::compute_hrr_hg_out_of_first(buffer, coordinates, 302578, 287878, 294178, 15,
                                             nmax);

        simdtrf::compute_hrr_hg_out_of_first(buffer, coordinates, 307303, 291028, 298378, 15,
                                             nmax);

        simdtrf::transform_g_inner(buffer, 312028, 307303, 21, 15, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 312028, 135, nmax);

        simdtrf::transform_g_inner(buffer, 312028, 302578, 21, 15, nmax);

        simdtrf::transform_h_outer(values + 1485 * nvalues + n * npairs, nvalues, buffer, 312028,
                                   135, nmax);
    }

    for (size_t m = 0; m < 2970; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
