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


#include "SimdThreeCenterElectronRepulsionRsRecISK.hpp"

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
#include "SimdTransformI.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_isk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_isk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 79724, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 390 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 79724, 77288, 2016, dimensions);

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
                                                            4, 5, 6, 7, 8, 9, 10, 11, 12, 13},
                                                            ncols, fj, i * nprim_b + j, fq,
                                                            omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 20, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13}, ncols,
                                                        fj, i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 34, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 37, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 40, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 43, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 46, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 49, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 52, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 55, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 58, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 61, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 64, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 67, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 70, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 73, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 76, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 79, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 82, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 85, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 88, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 91, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 94, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 97, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 100, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 103, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 106, 0, 3, 7, 8,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 112, 0, 3, 8, 9,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 118, 0, 3, 9, 10,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 124, 0, 3, 10, 11,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 130, 0, 3, 11, 12,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 136, 0, 3, 12, 13,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 142, 0, 3, 13, 14,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 148, 0, 3, 14, 15,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 154, 0, 3, 15, 16,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 160, 0, 3, 16, 17,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 166, 0, 3, 17, 18,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 172, 0, 3, 21, 22,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 178, 0, 3, 22, 23,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 184, 0, 3, 23, 24,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 190, 0, 3, 24, 25,
                                                                       79, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 196, 0, 3, 25, 26,
                                                                       82, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 202, 0, 3, 26, 27,
                                                                       85, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 208, 0, 3, 27, 28,
                                                                       88, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 214, 0, 3, 28, 29,
                                                                       91, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 220, 0, 3, 29, 30,
                                                                       94, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 226, 0, 3, 30, 31,
                                                                       97, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 232, 0, 3, 31, 32,
                                                                       100, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 238, 0, 3, 34, 37,
                                                                       106, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 248, 0, 3, 37, 40,
                                                                       112, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 258, 0, 3, 40, 43,
                                                                       118, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 268, 0, 3, 43, 46,
                                                                       124, 130, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 278, 0, 3, 46, 49,
                                                                       130, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 288, 0, 3, 49, 52,
                                                                       136, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 298, 0, 3, 52, 55,
                                                                       142, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 308, 0, 3, 55, 58,
                                                                       148, 154, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 318, 0, 3, 58, 61,
                                                                       154, 160, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 328, 0, 3, 61, 64,
                                                                       160, 166, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 338, 0, 3, 70, 73,
                                                                       172, 178, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 348, 0, 3, 73, 76,
                                                                       178, 184, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 76, 79,
                                                                       184, 190, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 368, 0, 3, 79, 82,
                                                                       190, 196, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 378, 0, 3, 82, 85,
                                                                       196, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 388, 0, 3, 85, 88,
                                                                       202, 208, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 398, 0, 3, 88, 91,
                                                                       208, 214, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 408, 0, 3, 91, 94,
                                                                       214, 220, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 418, 0, 3, 94, 97,
                                                                       220, 226, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 97,
                                                                       100, 226, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 438, 0, 3, 106,
                                                                       112, 238, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 453, 0, 3, 112,
                                                                       118, 248, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 468, 0, 3, 118,
                                                                       124, 258, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 483, 0, 3, 124,
                                                                       130, 268, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 498, 0, 3, 130,
                                                                       136, 278, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 513, 0, 3, 136,
                                                                       142, 288, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 528, 0, 3, 142,
                                                                       148, 298, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 543, 0, 3, 148,
                                                                       154, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 558, 0, 3, 154,
                                                                       160, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 573, 0, 3, 172,
                                                                       178, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 588, 0, 3, 178,
                                                                       184, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 603, 0, 3, 184,
                                                                       190, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 618, 0, 3, 190,
                                                                       196, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 633, 0, 3, 196,
                                                                       202, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 648, 0, 3, 202,
                                                                       208, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 663, 0, 3, 208,
                                                                       214, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 678, 0, 3, 214,
                                                                       220, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 693, 0, 3, 220,
                                                                       226, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 708, 0, 3, 238,
                                                                       248, 438, 453, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 729, 0, 3, 248,
                                                                       258, 453, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 750, 0, 3, 258,
                                                                       268, 468, 483, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 771, 0, 3, 268,
                                                                       278, 483, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 792, 0, 3, 278,
                                                                       288, 498, 513, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 813, 0, 3, 288,
                                                                       298, 513, 528, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 834, 0, 3, 298,
                                                                       308, 528, 543, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 855, 0, 3, 308,
                                                                       318, 543, 558, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 876, 0, 3, 338,
                                                                       348, 573, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 897, 0, 3, 348,
                                                                       358, 588, 603, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 918, 0, 3, 358,
                                                                       368, 603, 618, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 939, 0, 3, 368,
                                                                       378, 618, 633, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 960, 0, 3, 378,
                                                                       388, 633, 648, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 981, 0, 3, 388,
                                                                       398, 648, 663, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1002, 0, 3, 398,
                                                                       408, 663, 678, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1023, 0, 3, 408,
                                                                       418, 678, 693, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1044, 0, 3, 438,
                                                                       453, 708, 729, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1072, 0, 3, 453,
                                                                       468, 729, 750, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1100, 0, 3, 468,
                                                                       483, 750, 771, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1128, 0, 3, 483,
                                                                       498, 771, 792, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1156, 0, 3, 498,
                                                                       513, 792, 813, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1184, 0, 3, 513,
                                                                       528, 813, 834, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1212, 0, 3, 528,
                                                                       543, 834, 855, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1240, 0, 3, 573,
                                                                       588, 876, 897, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 588,
                                                                       603, 897, 918, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1296, 0, 3, 603,
                                                                       618, 918, 939, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1324, 0, 3, 618,
                                                                       633, 939, 960, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1352, 0, 3, 633,
                                                                       648, 960, 981, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 648,
                                                                       663, 981, 1002, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1408, 0, 3, 663,
                                                                       678, 1002, 1023, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1436, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1439, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1442, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1445, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1448, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1451, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1454, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1457, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1460, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1463, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1466, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1469, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1472, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1475, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1478, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1481, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1484, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1487, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1490, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1493, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1496, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1499, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1502, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1505, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1508, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1511, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1514, 3, 7, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1523, 3, 8, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1532, 3, 9, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1541, 3, 10, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1550, 3, 11, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1559, 3, 12, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1568, 3, 13, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1577, 3, 14, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1586, 3, 15, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1595, 3, 16, 61,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1604, 3, 17, 64,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1613, 3, 18, 67,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1622, 3, 21, 70,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1631, 3, 22, 73,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1640, 3, 23, 76,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1649, 3, 24, 79,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1658, 3, 25, 82,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1667, 3, 26, 85,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1676, 3, 27, 88,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1685, 3, 28, 91,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1694, 3, 29, 94,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1703, 3, 30, 97,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1712, 3, 31, 100,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1721, 3, 32, 103,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1730, 3, 34, 106,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1748, 3, 37, 112,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1766, 3, 40, 118,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1784, 3, 43, 124,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1802, 3, 46, 130,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1820, 3, 49, 136,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1838, 3, 52, 142,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1856, 3, 55, 148,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1874, 3, 58, 154,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1892, 3, 61, 160,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1910, 3, 64, 166,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1928, 3, 70, 172,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1946, 3, 73, 178,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1964, 3, 76, 184,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1982, 3, 79, 190,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2000, 3, 82, 196,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2018, 3, 85, 202,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2036, 3, 88, 208,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2054, 3, 91, 214,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2072, 3, 94, 220,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2090, 3, 97, 226,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2108, 3, 100, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2126, 3, 106, 238,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2156, 3, 112, 248,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2186, 3, 118, 258,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2216, 3, 124, 268,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2246, 3, 130, 278,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2276, 3, 136, 288,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2306, 3, 142, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2336, 3, 148, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2366, 3, 154, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2396, 3, 160, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2426, 3, 172, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2456, 3, 178, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2486, 3, 184, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2516, 3, 190, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2546, 3, 196, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2576, 3, 202, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2606, 3, 208, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2636, 3, 214, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2666, 3, 220, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2696, 3, 226, 428,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2726, 3, 238, 438,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2771, 3, 248, 453,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2816, 3, 258, 468,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2861, 3, 268, 483,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2906, 3, 278, 498,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2951, 3, 288, 513,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2996, 3, 298, 528,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3041, 3, 308, 543,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3086, 3, 318, 558,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3131, 3, 338, 573,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3176, 3, 348, 588,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3221, 3, 358, 603,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3266, 3, 368, 618,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3311, 3, 378, 633,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3356, 3, 388, 648,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3401, 3, 398, 663,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3446, 3, 408, 678,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3491, 3, 418, 693,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3536, 3, 438, 708,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3599, 3, 453, 729,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3662, 3, 468, 750,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3725, 3, 483, 771,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3788, 3, 498, 792,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3851, 3, 513, 813,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3914, 3, 528, 834,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3977, 3, 543, 855,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4040, 3, 573, 876,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4103, 3, 588, 897,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4166, 3, 603, 918,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4229, 3, 618, 939,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4292, 3, 633, 960,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4355, 3, 648, 981,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4418, 3, 663,
                                                                       1002, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4481, 3, 678,
                                                                       1023, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4544, 3, 708,
                                                                       1044, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4628, 3, 729,
                                                                       1072, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4712, 3, 750,
                                                                       1100, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4796, 3, 771,
                                                                       1128, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4880, 3, 792,
                                                                       1156, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4964, 3, 813,
                                                                       1184, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5048, 3, 834,
                                                                       1212, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5132, 3, 876,
                                                                       1240, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5216, 3, 897,
                                                                       1268, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5300, 3, 918,
                                                                       1296, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5384, 3, 939,
                                                                       1324, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5468, 3, 960,
                                                                       1352, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5552, 3, 981,
                                                                       1380, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5636, 3, 1002,
                                                                       1408, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5720, 3, 7, 8,
                                                                       1442, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5726, 3, 8, 9,
                                                                       1445, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5732, 3, 9, 10,
                                                                       1448, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5738, 3, 10, 11,
                                                                       1451, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5744, 3, 11, 12,
                                                                       1454, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5750, 3, 12, 13,
                                                                       1457, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5756, 3, 13, 14,
                                                                       1460, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5762, 3, 14, 15,
                                                                       1463, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5768, 3, 15, 16,
                                                                       1466, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5774, 3, 16, 17,
                                                                       1469, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5780, 3, 17, 18,
                                                                       1472, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5786, 3, 21, 22,
                                                                       1481, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5792, 3, 22, 23,
                                                                       1484, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5798, 3, 23, 24,
                                                                       1487, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5804, 3, 24, 25,
                                                                       1490, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5810, 3, 25, 26,
                                                                       1493, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5816, 3, 26, 27,
                                                                       1496, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5822, 3, 27, 28,
                                                                       1499, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5828, 3, 28, 29,
                                                                       1502, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5834, 3, 29, 30,
                                                                       1505, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5840, 3, 30, 31,
                                                                       1508, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5846, 3, 31, 32,
                                                                       1511, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5852, 0, 3, 5720,
                                                                       1442, 5726, 1532, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5870, 0, 3, 5726,
                                                                       1445, 5732, 1541, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5888, 0, 3, 5732,
                                                                       1448, 5738, 1550, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5906, 0, 3, 5738,
                                                                       1451, 5744, 1559, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5924, 0, 3, 5744,
                                                                       1454, 5750, 1568, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5942, 0, 3, 5750,
                                                                       1457, 5756, 1577, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5960, 0, 3, 5756,
                                                                       1460, 5762, 1586, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5978, 0, 3, 5762,
                                                                       1463, 5768, 1595, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5996, 0, 3, 5768,
                                                                       1466, 5774, 1604, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6014, 0, 3, 5774,
                                                                       1469, 5780, 1613, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6032, 0, 3, 5786,
                                                                       1481, 5792, 1640, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6050, 0, 3, 5792,
                                                                       1484, 5798, 1649, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6068, 0, 3, 5798,
                                                                       1487, 5804, 1658, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6086, 0, 3, 5804,
                                                                       1490, 5810, 1667, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6104, 0, 3, 5810,
                                                                       1493, 5816, 1676, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6122, 0, 3, 5816,
                                                                       1496, 5822, 1685, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6140, 0, 3, 5822,
                                                                       1499, 5828, 1694, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6158, 0, 3, 5828,
                                                                       1502, 5834, 1703, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6176, 0, 3, 5834,
                                                                       1505, 5840, 1712, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6194, 0, 3, 5840,
                                                                       1508, 5846, 1721, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6212, 0, 3, 5852,
                                                                       1532, 5870, 106, 112,
                                                                       1766, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6248, 0, 3, 5870,
                                                                       1541, 5888, 112, 118,
                                                                       1784, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6284, 0, 3, 5888,
                                                                       1550, 5906, 118, 124,
                                                                       1802, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6320, 0, 3, 5906,
                                                                       1559, 5924, 124, 130,
                                                                       1820, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6356, 0, 3, 5924,
                                                                       1568, 5942, 130, 136,
                                                                       1838, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6392, 0, 3, 5942,
                                                                       1577, 5960, 136, 142,
                                                                       1856, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6428, 0, 3, 5960,
                                                                       1586, 5978, 142, 148,
                                                                       1874, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6464, 0, 3, 5978,
                                                                       1595, 5996, 148, 154,
                                                                       1892, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6500, 0, 3, 5996,
                                                                       1604, 6014, 154, 160,
                                                                       1910, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6536, 0, 3, 6032,
                                                                       1640, 6050, 172, 178,
                                                                       1964, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6572, 0, 3, 6050,
                                                                       1649, 6068, 178, 184,
                                                                       1982, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6608, 0, 3, 6068,
                                                                       1658, 6086, 184, 190,
                                                                       2000, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6644, 0, 3, 6086,
                                                                       1667, 6104, 190, 196,
                                                                       2018, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6680, 0, 3, 6104,
                                                                       1676, 6122, 196, 202,
                                                                       2036, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6716, 0, 3, 6122,
                                                                       1685, 6140, 202, 208,
                                                                       2054, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6752, 0, 3, 6140,
                                                                       1694, 6158, 208, 214,
                                                                       2072, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6788, 0, 3, 6158,
                                                                       1703, 6176, 214, 220,
                                                                       2090, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6824, 0, 3, 6176,
                                                                       1712, 6194, 220, 226,
                                                                       2108, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6860, 0, 3, 6212,
                                                                       1766, 6248, 238, 248,
                                                                       2186, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6920, 0, 3, 6248,
                                                                       1784, 6284, 248, 258,
                                                                       2216, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6980, 0, 3, 6284,
                                                                       1802, 6320, 258, 268,
                                                                       2246, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7040, 0, 3, 6320,
                                                                       1820, 6356, 268, 278,
                                                                       2276, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7100, 0, 3, 6356,
                                                                       1838, 6392, 278, 288,
                                                                       2306, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7160, 0, 3, 6392,
                                                                       1856, 6428, 288, 298,
                                                                       2336, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7220, 0, 3, 6428,
                                                                       1874, 6464, 298, 308,
                                                                       2366, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7280, 0, 3, 6464,
                                                                       1892, 6500, 308, 318,
                                                                       2396, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7340, 0, 3, 6536,
                                                                       1964, 6572, 338, 348,
                                                                       2486, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7400, 0, 3, 6572,
                                                                       1982, 6608, 348, 358,
                                                                       2516, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7460, 0, 3, 6608,
                                                                       2000, 6644, 358, 368,
                                                                       2546, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7520, 0, 3, 6644,
                                                                       2018, 6680, 368, 378,
                                                                       2576, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7580, 0, 3, 6680,
                                                                       2036, 6716, 378, 388,
                                                                       2606, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7640, 0, 3, 6716,
                                                                       2054, 6752, 388, 398,
                                                                       2636, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7700, 0, 3, 6752,
                                                                       2072, 6788, 398, 408,
                                                                       2666, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7760, 0, 3, 6788,
                                                                       2090, 6824, 408, 418,
                                                                       2696, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7820, 0, 3, 6860,
                                                                       2186, 6920, 438, 453,
                                                                       2816, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7910, 0, 3, 6920,
                                                                       2216, 6980, 453, 468,
                                                                       2861, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8000, 0, 3, 6980,
                                                                       2246, 7040, 468, 483,
                                                                       2906, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8090, 0, 3, 7040,
                                                                       2276, 7100, 483, 498,
                                                                       2951, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8180, 0, 3, 7100,
                                                                       2306, 7160, 498, 513,
                                                                       2996, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8270, 0, 3, 7160,
                                                                       2336, 7220, 513, 528,
                                                                       3041, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8360, 0, 3, 7220,
                                                                       2366, 7280, 528, 543,
                                                                       3086, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8450, 0, 3, 7340,
                                                                       2486, 7400, 573, 588,
                                                                       3221, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8540, 0, 3, 7400,
                                                                       2516, 7460, 588, 603,
                                                                       3266, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8630, 0, 3, 7460,
                                                                       2546, 7520, 603, 618,
                                                                       3311, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8720, 0, 3, 7520,
                                                                       2576, 7580, 618, 633,
                                                                       3356, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8810, 0, 3, 7580,
                                                                       2606, 7640, 633, 648,
                                                                       3401, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8900, 0, 3, 7640,
                                                                       2636, 7700, 648, 663,
                                                                       3446, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8990, 0, 3, 7700,
                                                                       2666, 7760, 663, 678,
                                                                       3491, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9080, 0, 3, 7820,
                                                                       2816, 7910, 708, 729,
                                                                       3662, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9206, 0, 3, 7910,
                                                                       2861, 8000, 729, 750,
                                                                       3725, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9332, 0, 3, 8000,
                                                                       2906, 8090, 750, 771,
                                                                       3788, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9458, 0, 3, 8090,
                                                                       2951, 8180, 771, 792,
                                                                       3851, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9584, 0, 3, 8180,
                                                                       2996, 8270, 792, 813,
                                                                       3914, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9710, 0, 3, 8270,
                                                                       3041, 8360, 813, 834,
                                                                       3977, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9836, 0, 3, 8450,
                                                                       3221, 8540, 876, 897,
                                                                       4166, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9962, 0, 3, 8540,
                                                                       3266, 8630, 897, 918,
                                                                       4229, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10088, 0, 3, 8630,
                                                                       3311, 8720, 918, 939,
                                                                       4292, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10214, 0, 3, 8720,
                                                                       3356, 8810, 939, 960,
                                                                       4355, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10340, 0, 3, 8810,
                                                                       3401, 8900, 960, 981,
                                                                       4418, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10466, 0, 3, 8900,
                                                                       3446, 8990, 981, 1002,
                                                                       4481, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10592, 0, 3, 9080,
                                                                       3662, 9206, 1044, 1072,
                                                                       4712, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10760, 0, 3, 9206,
                                                                       3725, 9332, 1072, 1100,
                                                                       4796, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10928, 0, 3, 9332,
                                                                       3788, 9458, 1100, 1128,
                                                                       4880, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11096, 0, 3, 9458,
                                                                       3851, 9584, 1128, 1156,
                                                                       4964, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11264, 0, 3, 9584,
                                                                       3914, 9710, 1156, 1184,
                                                                       5048, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11432, 0, 3, 9836,
                                                                       4166, 9962, 1240, 1268,
                                                                       5300, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11600, 0, 3, 9962,
                                                                       4229, 10088, 1268, 1296,
                                                                       5384, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11768, 0, 3,
                                                                       10088, 4292, 10214, 1296,
                                                                       1324, 5468, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11936, 0, 3,
                                                                       10214, 4355, 10340, 1324,
                                                                       1352, 5552, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12104, 0, 3,
                                                                       10340, 4418, 10466, 1352,
                                                                       1380, 5636, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12272, 3, 1436,
                                                                       1439, 5720, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12282, 3, 1439,
                                                                       1442, 5726, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12292, 3, 1442,
                                                                       1445, 5732, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12302, 3, 1445,
                                                                       1448, 5738, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12312, 3, 1448,
                                                                       1451, 5744, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12322, 3, 1451,
                                                                       1454, 5750, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12332, 3, 1454,
                                                                       1457, 5756, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12342, 3, 1457,
                                                                       1460, 5762, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12352, 3, 1460,
                                                                       1463, 5768, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12362, 3, 1463,
                                                                       1466, 5774, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12372, 3, 1466,
                                                                       1469, 5780, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12382, 3, 1475,
                                                                       1478, 5786, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12392, 3, 1478,
                                                                       1481, 5792, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12402, 3, 1481,
                                                                       1484, 5798, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12412, 3, 1484,
                                                                       1487, 5804, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12422, 3, 1487,
                                                                       1490, 5810, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12432, 3, 1490,
                                                                       1493, 5816, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12442, 3, 1493,
                                                                       1496, 5822, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12452, 3, 1496,
                                                                       1499, 5828, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12462, 3, 1499,
                                                                       1502, 5834, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12472, 3, 1502,
                                                                       1505, 5840, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12482, 3, 1505,
                                                                       1508, 5846, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12492, 0, 3,
                                                                       12272, 5720, 12282, 1514,
                                                                       1523, 5852, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12522, 0, 3,
                                                                       12282, 5726, 12292, 1523,
                                                                       1532, 5870, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12552, 0, 3,
                                                                       12292, 5732, 12302, 1532,
                                                                       1541, 5888, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12582, 0, 3,
                                                                       12302, 5738, 12312, 1541,
                                                                       1550, 5906, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12612, 0, 3,
                                                                       12312, 5744, 12322, 1550,
                                                                       1559, 5924, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12642, 0, 3,
                                                                       12322, 5750, 12332, 1559,
                                                                       1568, 5942, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12672, 0, 3,
                                                                       12332, 5756, 12342, 1568,
                                                                       1577, 5960, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12702, 0, 3,
                                                                       12342, 5762, 12352, 1577,
                                                                       1586, 5978, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12732, 0, 3,
                                                                       12352, 5768, 12362, 1586,
                                                                       1595, 5996, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12762, 0, 3,
                                                                       12362, 5774, 12372, 1595,
                                                                       1604, 6014, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12792, 0, 3,
                                                                       12382, 5786, 12392, 1622,
                                                                       1631, 6032, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12822, 0, 3,
                                                                       12392, 5792, 12402, 1631,
                                                                       1640, 6050, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12852, 0, 3,
                                                                       12402, 5798, 12412, 1640,
                                                                       1649, 6068, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12882, 0, 3,
                                                                       12412, 5804, 12422, 1649,
                                                                       1658, 6086, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12912, 0, 3,
                                                                       12422, 5810, 12432, 1658,
                                                                       1667, 6104, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12942, 0, 3,
                                                                       12432, 5816, 12442, 1667,
                                                                       1676, 6122, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12972, 0, 3,
                                                                       12442, 5822, 12452, 1676,
                                                                       1685, 6140, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13002, 0, 3,
                                                                       12452, 5828, 12462, 1685,
                                                                       1694, 6158, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13032, 0, 3,
                                                                       12462, 5834, 12472, 1694,
                                                                       1703, 6176, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13062, 0, 3,
                                                                       12472, 5840, 12482, 1703,
                                                                       1712, 6194, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13092, 0, 3,
                                                                       12492, 5852, 12522, 1730,
                                                                       1748, 6212, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13152, 0, 3,
                                                                       12522, 5870, 12552, 1748,
                                                                       1766, 6248, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13212, 0, 3,
                                                                       12552, 5888, 12582, 1766,
                                                                       1784, 6284, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13272, 0, 3,
                                                                       12582, 5906, 12612, 1784,
                                                                       1802, 6320, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13332, 0, 3,
                                                                       12612, 5924, 12642, 1802,
                                                                       1820, 6356, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13392, 0, 3,
                                                                       12642, 5942, 12672, 1820,
                                                                       1838, 6392, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13452, 0, 3,
                                                                       12672, 5960, 12702, 1838,
                                                                       1856, 6428, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13512, 0, 3,
                                                                       12702, 5978, 12732, 1856,
                                                                       1874, 6464, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13572, 0, 3,
                                                                       12732, 5996, 12762, 1874,
                                                                       1892, 6500, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13632, 0, 3,
                                                                       12792, 6032, 12822, 1928,
                                                                       1946, 6536, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13692, 0, 3,
                                                                       12822, 6050, 12852, 1946,
                                                                       1964, 6572, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13752, 0, 3,
                                                                       12852, 6068, 12882, 1964,
                                                                       1982, 6608, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13812, 0, 3,
                                                                       12882, 6086, 12912, 1982,
                                                                       2000, 6644, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13872, 0, 3,
                                                                       12912, 6104, 12942, 2000,
                                                                       2018, 6680, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13932, 0, 3,
                                                                       12942, 6122, 12972, 2018,
                                                                       2036, 6716, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13992, 0, 3,
                                                                       12972, 6140, 13002, 2036,
                                                                       2054, 6752, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14052, 0, 3,
                                                                       13002, 6158, 13032, 2054,
                                                                       2072, 6788, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14112, 0, 3,
                                                                       13032, 6176, 13062, 2072,
                                                                       2090, 6824, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14172, 0, 3,
                                                                       13092, 6212, 13152, 2126,
                                                                       2156, 6860, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14272, 0, 3,
                                                                       13152, 6248, 13212, 2156,
                                                                       2186, 6920, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14372, 0, 3,
                                                                       13212, 6284, 13272, 2186,
                                                                       2216, 6980, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14472, 0, 3,
                                                                       13272, 6320, 13332, 2216,
                                                                       2246, 7040, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14572, 0, 3,
                                                                       13332, 6356, 13392, 2246,
                                                                       2276, 7100, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14672, 0, 3,
                                                                       13392, 6392, 13452, 2276,
                                                                       2306, 7160, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14772, 0, 3,
                                                                       13452, 6428, 13512, 2306,
                                                                       2336, 7220, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14872, 0, 3,
                                                                       13512, 6464, 13572, 2336,
                                                                       2366, 7280, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14972, 0, 3,
                                                                       13632, 6536, 13692, 2426,
                                                                       2456, 7340, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15072, 0, 3,
                                                                       13692, 6572, 13752, 2456,
                                                                       2486, 7400, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15172, 0, 3,
                                                                       13752, 6608, 13812, 2486,
                                                                       2516, 7460, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15272, 0, 3,
                                                                       13812, 6644, 13872, 2516,
                                                                       2546, 7520, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15372, 0, 3,
                                                                       13872, 6680, 13932, 2546,
                                                                       2576, 7580, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15472, 0, 3,
                                                                       13932, 6716, 13992, 2576,
                                                                       2606, 7640, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15572, 0, 3,
                                                                       13992, 6752, 14052, 2606,
                                                                       2636, 7700, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15672, 0, 3,
                                                                       14052, 6788, 14112, 2636,
                                                                       2666, 7760, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 15772, 0, 3,
                                                                       14172, 6860, 14272, 2726,
                                                                       2771, 7820, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 15922, 0, 3,
                                                                       14272, 6920, 14372, 2771,
                                                                       2816, 7910, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 16072, 0, 3,
                                                                       14372, 6980, 14472, 2816,
                                                                       2861, 8000, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 16222, 0, 3,
                                                                       14472, 7040, 14572, 2861,
                                                                       2906, 8090, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 16372, 0, 3,
                                                                       14572, 7100, 14672, 2906,
                                                                       2951, 8180, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 16522, 0, 3,
                                                                       14672, 7160, 14772, 2951,
                                                                       2996, 8270, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 16672, 0, 3,
                                                                       14772, 7220, 14872, 2996,
                                                                       3041, 8360, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 16822, 0, 3,
                                                                       14972, 7340, 15072, 3131,
                                                                       3176, 8450, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 16972, 0, 3,
                                                                       15072, 7400, 15172, 3176,
                                                                       3221, 8540, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17122, 0, 3,
                                                                       15172, 7460, 15272, 3221,
                                                                       3266, 8630, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17272, 0, 3,
                                                                       15272, 7520, 15372, 3266,
                                                                       3311, 8720, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17422, 0, 3,
                                                                       15372, 7580, 15472, 3311,
                                                                       3356, 8810, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17572, 0, 3,
                                                                       15472, 7640, 15572, 3356,
                                                                       3401, 8900, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17722, 0, 3,
                                                                       15572, 7700, 15672, 3401,
                                                                       3446, 8990, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 17872, 0, 3,
                                                                       15772, 7820, 15922, 3536,
                                                                       3599, 9080, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 18082, 0, 3,
                                                                       15922, 7910, 16072, 3599,
                                                                       3662, 9206, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 18292, 0, 3,
                                                                       16072, 8000, 16222, 3662,
                                                                       3725, 9332, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 18502, 0, 3,
                                                                       16222, 8090, 16372, 3725,
                                                                       3788, 9458, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 18712, 0, 3,
                                                                       16372, 8180, 16522, 3788,
                                                                       3851, 9584, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 18922, 0, 3,
                                                                       16522, 8270, 16672, 3851,
                                                                       3914, 9710, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19132, 0, 3,
                                                                       16822, 8450, 16972, 4040,
                                                                       4103, 9836, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19342, 0, 3,
                                                                       16972, 8540, 17122, 4103,
                                                                       4166, 9962, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19552, 0, 3,
                                                                       17122, 8630, 17272, 4166,
                                                                       4229, 10088, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19762, 0, 3,
                                                                       17272, 8720, 17422, 4229,
                                                                       4292, 10214, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19972, 0, 3,
                                                                       17422, 8810, 17572, 4292,
                                                                       4355, 10340, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 20182, 0, 3,
                                                                       17572, 8900, 17722, 4355,
                                                                       4418, 10466, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 20392, 0, 3,
                                                                       17872, 9080, 18082, 4544,
                                                                       4628, 10592, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 20672, 0, 3,
                                                                       18082, 9206, 18292, 4628,
                                                                       4712, 10760, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 20952, 0, 3,
                                                                       18292, 9332, 18502, 4712,
                                                                       4796, 10928, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 21232, 0, 3,
                                                                       18502, 9458, 18712, 4796,
                                                                       4880, 11096, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 21512, 0, 3,
                                                                       18712, 9584, 18922, 4880,
                                                                       4964, 11264, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 21792, 0, 3,
                                                                       19132, 9836, 19342, 5132,
                                                                       5216, 11432, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 22072, 0, 3,
                                                                       19342, 9962, 19552, 5216,
                                                                       5300, 11600, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 22352, 0, 3,
                                                                       19552, 10088, 19762, 5300,
                                                                       5384, 11768, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 22632, 0, 3,
                                                                       19762, 10214, 19972, 5384,
                                                                       5468, 11936, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 22912, 0, 3,
                                                                       19972, 10340, 20182, 5468,
                                                                       5552, 12104, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23192, 3, 5720,
                                                                       5726, 12292, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23207, 3, 5726,
                                                                       5732, 12302, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23222, 3, 5732,
                                                                       5738, 12312, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23237, 3, 5738,
                                                                       5744, 12322, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23252, 3, 5744,
                                                                       5750, 12332, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23267, 3, 5750,
                                                                       5756, 12342, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23282, 3, 5756,
                                                                       5762, 12352, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23297, 3, 5762,
                                                                       5768, 12362, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23312, 3, 5768,
                                                                       5774, 12372, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23327, 3, 5786,
                                                                       5792, 12402, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23342, 3, 5792,
                                                                       5798, 12412, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23357, 3, 5798,
                                                                       5804, 12422, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23372, 3, 5804,
                                                                       5810, 12432, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23387, 3, 5810,
                                                                       5816, 12442, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23402, 3, 5816,
                                                                       5822, 12452, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23417, 3, 5822,
                                                                       5828, 12462, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23432, 3, 5828,
                                                                       5834, 12472, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23447, 3, 5834,
                                                                       5840, 12482, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 23462, 0, 3,
                                                                       23192, 12292, 23207, 5852,
                                                                       5870, 12552, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 23507, 0, 3,
                                                                       23207, 12302, 23222, 5870,
                                                                       5888, 12582, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 23552, 0, 3,
                                                                       23222, 12312, 23237, 5888,
                                                                       5906, 12612, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 23597, 0, 3,
                                                                       23237, 12322, 23252, 5906,
                                                                       5924, 12642, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 23642, 0, 3,
                                                                       23252, 12332, 23267, 5924,
                                                                       5942, 12672, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 23687, 0, 3,
                                                                       23267, 12342, 23282, 5942,
                                                                       5960, 12702, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 23732, 0, 3,
                                                                       23282, 12352, 23297, 5960,
                                                                       5978, 12732, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 23777, 0, 3,
                                                                       23297, 12362, 23312, 5978,
                                                                       5996, 12762, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 23822, 0, 3,
                                                                       23327, 12402, 23342, 6032,
                                                                       6050, 12852, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 23867, 0, 3,
                                                                       23342, 12412, 23357, 6050,
                                                                       6068, 12882, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 23912, 0, 3,
                                                                       23357, 12422, 23372, 6068,
                                                                       6086, 12912, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 23957, 0, 3,
                                                                       23372, 12432, 23387, 6086,
                                                                       6104, 12942, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24002, 0, 3,
                                                                       23387, 12442, 23402, 6104,
                                                                       6122, 12972, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24047, 0, 3,
                                                                       23402, 12452, 23417, 6122,
                                                                       6140, 13002, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24092, 0, 3,
                                                                       23417, 12462, 23432, 6140,
                                                                       6158, 13032, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24137, 0, 3,
                                                                       23432, 12472, 23447, 6158,
                                                                       6176, 13062, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 24182, 0, 3,
                                                                       23462, 12552, 23507, 6212,
                                                                       6248, 13212, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 24272, 0, 3,
                                                                       23507, 12582, 23552, 6248,
                                                                       6284, 13272, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 24362, 0, 3,
                                                                       23552, 12612, 23597, 6284,
                                                                       6320, 13332, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 24452, 0, 3,
                                                                       23597, 12642, 23642, 6320,
                                                                       6356, 13392, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 24542, 0, 3,
                                                                       23642, 12672, 23687, 6356,
                                                                       6392, 13452, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 24632, 0, 3,
                                                                       23687, 12702, 23732, 6392,
                                                                       6428, 13512, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 24722, 0, 3,
                                                                       23732, 12732, 23777, 6428,
                                                                       6464, 13572, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 24812, 0, 3,
                                                                       23822, 12852, 23867, 6536,
                                                                       6572, 13752, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 24902, 0, 3,
                                                                       23867, 12882, 23912, 6572,
                                                                       6608, 13812, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 24992, 0, 3,
                                                                       23912, 12912, 23957, 6608,
                                                                       6644, 13872, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25082, 0, 3,
                                                                       23957, 12942, 24002, 6644,
                                                                       6680, 13932, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25172, 0, 3,
                                                                       24002, 12972, 24047, 6680,
                                                                       6716, 13992, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25262, 0, 3,
                                                                       24047, 13002, 24092, 6716,
                                                                       6752, 14052, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25352, 0, 3,
                                                                       24092, 13032, 24137, 6752,
                                                                       6788, 14112, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 25442, 0, 3,
                                                                       24182, 13212, 24272, 6860,
                                                                       6920, 14372, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 25592, 0, 3,
                                                                       24272, 13272, 24362, 6920,
                                                                       6980, 14472, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 25742, 0, 3,
                                                                       24362, 13332, 24452, 6980,
                                                                       7040, 14572, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 25892, 0, 3,
                                                                       24452, 13392, 24542, 7040,
                                                                       7100, 14672, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26042, 0, 3,
                                                                       24542, 13452, 24632, 7100,
                                                                       7160, 14772, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26192, 0, 3,
                                                                       24632, 13512, 24722, 7160,
                                                                       7220, 14872, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26342, 0, 3,
                                                                       24812, 13752, 24902, 7340,
                                                                       7400, 15172, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26492, 0, 3,
                                                                       24902, 13812, 24992, 7400,
                                                                       7460, 15272, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26642, 0, 3,
                                                                       24992, 13872, 25082, 7460,
                                                                       7520, 15372, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26792, 0, 3,
                                                                       25082, 13932, 25172, 7520,
                                                                       7580, 15472, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26942, 0, 3,
                                                                       25172, 13992, 25262, 7580,
                                                                       7640, 15572, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 27092, 0, 3,
                                                                       25262, 14052, 25352, 7640,
                                                                       7700, 15672, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 27242, 0, 3,
                                                                       25442, 14372, 25592, 7820,
                                                                       7910, 16072, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 27467, 0, 3,
                                                                       25592, 14472, 25742, 7910,
                                                                       8000, 16222, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 27692, 0, 3,
                                                                       25742, 14572, 25892, 8000,
                                                                       8090, 16372, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 27917, 0, 3,
                                                                       25892, 14672, 26042, 8090,
                                                                       8180, 16522, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 28142, 0, 3,
                                                                       26042, 14772, 26192, 8180,
                                                                       8270, 16672, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 28367, 0, 3,
                                                                       26342, 15172, 26492, 8450,
                                                                       8540, 17122, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 28592, 0, 3,
                                                                       26492, 15272, 26642, 8540,
                                                                       8630, 17272, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 28817, 0, 3,
                                                                       26642, 15372, 26792, 8630,
                                                                       8720, 17422, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 29042, 0, 3,
                                                                       26792, 15472, 26942, 8720,
                                                                       8810, 17572, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 29267, 0, 3,
                                                                       26942, 15572, 27092, 8810,
                                                                       8900, 17722, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 29492, 0, 3,
                                                                       27242, 16072, 27467, 9080,
                                                                       9206, 18292, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 29807, 0, 3,
                                                                       27467, 16222, 27692, 9206,
                                                                       9332, 18502, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 30122, 0, 3,
                                                                       27692, 16372, 27917, 9332,
                                                                       9458, 18712, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 30437, 0, 3,
                                                                       27917, 16522, 28142, 9458,
                                                                       9584, 18922, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 30752, 0, 3,
                                                                       28367, 17122, 28592, 9836,
                                                                       9962, 19552, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 31067, 0, 3,
                                                                       28592, 17272, 28817, 9962,
                                                                       10088, 19762, ncols,
                                                                       gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 31382, 0, 3,
                                                                       28817, 17422, 29042,
                                                                       10088, 10214, 19972,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 31697, 0, 3,
                                                                       29042, 17572, 29267,
                                                                       10214, 10340, 20182,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 32012, 0, 3,
                                                                       29492, 18292, 29807,
                                                                       10592, 10760, 20952,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 32432, 0, 3,
                                                                       29807, 18502, 30122,
                                                                       10760, 10928, 21232,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 32852, 0, 3,
                                                                       30122, 18712, 30437,
                                                                       10928, 11096, 21512,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 33272, 0, 3,
                                                                       30752, 19552, 31067,
                                                                       11432, 11600, 22352,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 33692, 0, 3,
                                                                       31067, 19762, 31382,
                                                                       11600, 11768, 22632,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 34112, 0, 3,
                                                                       31382, 19972, 31697,
                                                                       11768, 11936, 22912,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34532, 3, 12272,
                                                                       12282, 23192, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34553, 3, 12282,
                                                                       12292, 23207, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34574, 3, 12292,
                                                                       12302, 23222, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34595, 3, 12302,
                                                                       12312, 23237, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34616, 3, 12312,
                                                                       12322, 23252, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34637, 3, 12322,
                                                                       12332, 23267, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34658, 3, 12332,
                                                                       12342, 23282, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34679, 3, 12342,
                                                                       12352, 23297, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34700, 3, 12352,
                                                                       12362, 23312, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34721, 3, 12382,
                                                                       12392, 23327, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34742, 3, 12392,
                                                                       12402, 23342, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34763, 3, 12402,
                                                                       12412, 23357, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34784, 3, 12412,
                                                                       12422, 23372, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34805, 3, 12422,
                                                                       12432, 23387, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34826, 3, 12432,
                                                                       12442, 23402, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34847, 3, 12442,
                                                                       12452, 23417, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34868, 3, 12452,
                                                                       12462, 23432, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34889, 3, 12462,
                                                                       12472, 23447, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 34910, 0, 3,
                                                                       34532, 23192, 34553,
                                                                       12492, 12522, 23462,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 34973, 0, 3,
                                                                       34553, 23207, 34574,
                                                                       12522, 12552, 23507,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35036, 0, 3,
                                                                       34574, 23222, 34595,
                                                                       12552, 12582, 23552,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35099, 0, 3,
                                                                       34595, 23237, 34616,
                                                                       12582, 12612, 23597,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35162, 0, 3,
                                                                       34616, 23252, 34637,
                                                                       12612, 12642, 23642,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35225, 0, 3,
                                                                       34637, 23267, 34658,
                                                                       12642, 12672, 23687,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35288, 0, 3,
                                                                       34658, 23282, 34679,
                                                                       12672, 12702, 23732,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35351, 0, 3,
                                                                       34679, 23297, 34700,
                                                                       12702, 12732, 23777,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35414, 0, 3,
                                                                       34721, 23327, 34742,
                                                                       12792, 12822, 23822,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35477, 0, 3,
                                                                       34742, 23342, 34763,
                                                                       12822, 12852, 23867,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35540, 0, 3,
                                                                       34763, 23357, 34784,
                                                                       12852, 12882, 23912,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35603, 0, 3,
                                                                       34784, 23372, 34805,
                                                                       12882, 12912, 23957,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35666, 0, 3,
                                                                       34805, 23387, 34826,
                                                                       12912, 12942, 24002,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35729, 0, 3,
                                                                       34826, 23402, 34847,
                                                                       12942, 12972, 24047,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35792, 0, 3,
                                                                       34847, 23417, 34868,
                                                                       12972, 13002, 24092,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35855, 0, 3,
                                                                       34868, 23432, 34889,
                                                                       13002, 13032, 24137,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 35918, 0, 3,
                                                                       34910, 23462, 34973,
                                                                       13092, 13152, 24182,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 36044, 0, 3,
                                                                       34973, 23507, 35036,
                                                                       13152, 13212, 24272,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 36170, 0, 3,
                                                                       35036, 23552, 35099,
                                                                       13212, 13272, 24362,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 36296, 0, 3,
                                                                       35099, 23597, 35162,
                                                                       13272, 13332, 24452,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 36422, 0, 3,
                                                                       35162, 23642, 35225,
                                                                       13332, 13392, 24542,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 36548, 0, 3,
                                                                       35225, 23687, 35288,
                                                                       13392, 13452, 24632,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 36674, 0, 3,
                                                                       35288, 23732, 35351,
                                                                       13452, 13512, 24722,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 36800, 0, 3,
                                                                       35414, 23822, 35477,
                                                                       13632, 13692, 24812,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 36926, 0, 3,
                                                                       35477, 23867, 35540,
                                                                       13692, 13752, 24902,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 37052, 0, 3,
                                                                       35540, 23912, 35603,
                                                                       13752, 13812, 24992,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 37178, 0, 3,
                                                                       35603, 23957, 35666,
                                                                       13812, 13872, 25082,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 37304, 0, 3,
                                                                       35666, 24002, 35729,
                                                                       13872, 13932, 25172,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 37430, 0, 3,
                                                                       35729, 24047, 35792,
                                                                       13932, 13992, 25262,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 37556, 0, 3,
                                                                       35792, 24092, 35855,
                                                                       13992, 14052, 25352,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 37682, 0, 3,
                                                                       35918, 24182, 36044,
                                                                       14172, 14272, 25442,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 37892, 0, 3,
                                                                       36044, 24272, 36170,
                                                                       14272, 14372, 25592,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 38102, 0, 3,
                                                                       36170, 24362, 36296,
                                                                       14372, 14472, 25742,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 38312, 0, 3,
                                                                       36296, 24452, 36422,
                                                                       14472, 14572, 25892,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 38522, 0, 3,
                                                                       36422, 24542, 36548,
                                                                       14572, 14672, 26042,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 38732, 0, 3,
                                                                       36548, 24632, 36674,
                                                                       14672, 14772, 26192,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 38942, 0, 3,
                                                                       36800, 24812, 36926,
                                                                       14972, 15072, 26342,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 39152, 0, 3,
                                                                       36926, 24902, 37052,
                                                                       15072, 15172, 26492,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 39362, 0, 3,
                                                                       37052, 24992, 37178,
                                                                       15172, 15272, 26642,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 39572, 0, 3,
                                                                       37178, 25082, 37304,
                                                                       15272, 15372, 26792,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 39782, 0, 3,
                                                                       37304, 25172, 37430,
                                                                       15372, 15472, 26942,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 39992, 0, 3,
                                                                       37430, 25262, 37556,
                                                                       15472, 15572, 27092,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 40202, 0, 3,
                                                                       37682, 25442, 37892,
                                                                       15772, 15922, 27242,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 40517, 0, 3,
                                                                       37892, 25592, 38102,
                                                                       15922, 16072, 27467,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 40832, 0, 3,
                                                                       38102, 25742, 38312,
                                                                       16072, 16222, 27692,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 41147, 0, 3,
                                                                       38312, 25892, 38522,
                                                                       16222, 16372, 27917,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 41462, 0, 3,
                                                                       38522, 26042, 38732,
                                                                       16372, 16522, 28142,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 41777, 0, 3,
                                                                       38942, 26342, 39152,
                                                                       16822, 16972, 28367,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 42092, 0, 3,
                                                                       39152, 26492, 39362,
                                                                       16972, 17122, 28592,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 42407, 0, 3,
                                                                       39362, 26642, 39572,
                                                                       17122, 17272, 28817,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 42722, 0, 3,
                                                                       39572, 26792, 39782,
                                                                       17272, 17422, 29042,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 43037, 0, 3,
                                                                       39782, 26942, 39992,
                                                                       17422, 17572, 29267,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 43352, 0, 3,
                                                                       40202, 27242, 40517,
                                                                       17872, 18082, 29492,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 43793, 0, 3,
                                                                       40517, 27467, 40832,
                                                                       18082, 18292, 29807,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 44234, 0, 3,
                                                                       40832, 27692, 41147,
                                                                       18292, 18502, 30122,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 44675, 0, 3,
                                                                       41147, 27917, 41462,
                                                                       18502, 18712, 30437,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 45116, 0, 3,
                                                                       41777, 28367, 42092,
                                                                       19132, 19342, 30752,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 45557, 0, 3,
                                                                       42092, 28592, 42407,
                                                                       19342, 19552, 31067,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 45998, 0, 3,
                                                                       42407, 28817, 42722,
                                                                       19552, 19762, 31382,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 46439, 0, 3,
                                                                       42722, 29042, 43037,
                                                                       19762, 19972, 31697,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 46880, 0, 3,
                                                                       43352, 29492, 43793,
                                                                       20392, 20672, 32012,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 47468, 0, 3,
                                                                       43793, 29807, 44234,
                                                                       20672, 20952, 32432,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 48056, 0, 3,
                                                                       44234, 30122, 44675,
                                                                       20952, 21232, 32852,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 48644, 0, 3,
                                                                       45116, 30752, 45557,
                                                                       21792, 22072, 33272,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 49232, 0, 3,
                                                                       45557, 31067, 45998,
                                                                       22072, 22352, 33692,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 49820, 0, 3,
                                                                       45998, 31382, 46439,
                                                                       22352, 22632, 34112,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50408, 3, 23192,
                                                                       23207, 34574, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50436, 3, 23207,
                                                                       23222, 34595, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50464, 3, 23222,
                                                                       23237, 34616, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50492, 3, 23237,
                                                                       23252, 34637, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50520, 3, 23252,
                                                                       23267, 34658, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50548, 3, 23267,
                                                                       23282, 34679, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50576, 3, 23282,
                                                                       23297, 34700, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50604, 3, 23327,
                                                                       23342, 34763, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50632, 3, 23342,
                                                                       23357, 34784, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50660, 3, 23357,
                                                                       23372, 34805, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50688, 3, 23372,
                                                                       23387, 34826, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50716, 3, 23387,
                                                                       23402, 34847, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50744, 3, 23402,
                                                                       23417, 34868, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50772, 3, 23417,
                                                                       23432, 34889, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 50800, 0, 3,
                                                                       50408, 34574, 50436,
                                                                       23462, 23507, 35036,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 50884, 0, 3,
                                                                       50436, 34595, 50464,
                                                                       23507, 23552, 35099,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 50968, 0, 3,
                                                                       50464, 34616, 50492,
                                                                       23552, 23597, 35162,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 51052, 0, 3,
                                                                       50492, 34637, 50520,
                                                                       23597, 23642, 35225,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 51136, 0, 3,
                                                                       50520, 34658, 50548,
                                                                       23642, 23687, 35288,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 51220, 0, 3,
                                                                       50548, 34679, 50576,
                                                                       23687, 23732, 35351,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 51304, 0, 3,
                                                                       50604, 34763, 50632,
                                                                       23822, 23867, 35540,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 51388, 0, 3,
                                                                       50632, 34784, 50660,
                                                                       23867, 23912, 35603,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 51472, 0, 3,
                                                                       50660, 34805, 50688,
                                                                       23912, 23957, 35666,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 51556, 0, 3,
                                                                       50688, 34826, 50716,
                                                                       23957, 24002, 35729,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 51640, 0, 3,
                                                                       50716, 34847, 50744,
                                                                       24002, 24047, 35792,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 51724, 0, 3,
                                                                       50744, 34868, 50772,
                                                                       24047, 24092, 35855,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 51808, 0, 3,
                                                                       50800, 35036, 50884,
                                                                       24182, 24272, 36170,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 51976, 0, 3,
                                                                       50884, 35099, 50968,
                                                                       24272, 24362, 36296,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 52144, 0, 3,
                                                                       50968, 35162, 51052,
                                                                       24362, 24452, 36422,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 52312, 0, 3,
                                                                       51052, 35225, 51136,
                                                                       24452, 24542, 36548,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 52480, 0, 3,
                                                                       51136, 35288, 51220,
                                                                       24542, 24632, 36674,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 52648, 0, 3,
                                                                       51304, 35540, 51388,
                                                                       24812, 24902, 37052,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 52816, 0, 3,
                                                                       51388, 35603, 51472,
                                                                       24902, 24992, 37178,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 52984, 0, 3,
                                                                       51472, 35666, 51556,
                                                                       24992, 25082, 37304,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 53152, 0, 3,
                                                                       51556, 35729, 51640,
                                                                       25082, 25172, 37430,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 53320, 0, 3,
                                                                       51640, 35792, 51724,
                                                                       25172, 25262, 37556,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 53488, 0, 3,
                                                                       51808, 36170, 51976,
                                                                       25442, 25592, 38102,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 53768, 0, 3,
                                                                       51976, 36296, 52144,
                                                                       25592, 25742, 38312,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 54048, 0, 3,
                                                                       52144, 36422, 52312,
                                                                       25742, 25892, 38522,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 54328, 0, 3,
                                                                       52312, 36548, 52480,
                                                                       25892, 26042, 38732,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 54608, 0, 3,
                                                                       52648, 37052, 52816,
                                                                       26342, 26492, 39362,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 54888, 0, 3,
                                                                       52816, 37178, 52984,
                                                                       26492, 26642, 39572,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 55168, 0, 3,
                                                                       52984, 37304, 53152,
                                                                       26642, 26792, 39782,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 55448, 0, 3,
                                                                       53152, 37430, 53320,
                                                                       26792, 26942, 39992,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 55728, 0, 3,
                                                                       53488, 38102, 53768,
                                                                       27242, 27467, 40832,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 56148, 0, 3,
                                                                       53768, 38312, 54048,
                                                                       27467, 27692, 41147,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 56568, 0, 3,
                                                                       54048, 38522, 54328,
                                                                       27692, 27917, 41462,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 56988, 0, 3,
                                                                       54608, 39362, 54888,
                                                                       28367, 28592, 42407,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 57408, 0, 3,
                                                                       54888, 39572, 55168,
                                                                       28592, 28817, 42722,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 57828, 0, 3,
                                                                       55168, 39782, 55448,
                                                                       28817, 29042, 43037,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 58248, 0, 3,
                                                                       55728, 40832, 56148,
                                                                       29492, 29807, 44234,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 58836, 0, 3,
                                                                       56148, 41147, 56568,
                                                                       29807, 30122, 44675,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 59424, 0, 3,
                                                                       56988, 42407, 57408,
                                                                       30752, 31067, 45998,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 60012, 0, 3,
                                                                       57408, 42722, 57828,
                                                                       31067, 31382, 46439,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 60600, 0, 3,
                                                                       58248, 44234, 58836,
                                                                       32012, 32432, 48056,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 61384, 0, 3,
                                                                       59424, 45998, 60012,
                                                                       33272, 33692, 49820,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62168, 3, 34532,
                                                                       34553, 50408, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62204, 3, 34553,
                                                                       34574, 50436, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62240, 3, 34574,
                                                                       34595, 50464, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62276, 3, 34595,
                                                                       34616, 50492, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62312, 3, 34616,
                                                                       34637, 50520, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62348, 3, 34637,
                                                                       34658, 50548, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62384, 3, 34658,
                                                                       34679, 50576, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62420, 3, 34721,
                                                                       34742, 50604, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62456, 3, 34742,
                                                                       34763, 50632, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62492, 3, 34763,
                                                                       34784, 50660, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62528, 3, 34784,
                                                                       34805, 50688, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62564, 3, 34805,
                                                                       34826, 50716, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62600, 3, 34826,
                                                                       34847, 50744, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62636, 3, 34847,
                                                                       34868, 50772, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 62672, 0, 3,
                                                                       62168, 50408, 62204,
                                                                       34910, 34973, 50800,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 62780, 0, 3,
                                                                       62204, 50436, 62240,
                                                                       34973, 35036, 50884,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 62888, 0, 3,
                                                                       62240, 50464, 62276,
                                                                       35036, 35099, 50968,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 62996, 0, 3,
                                                                       62276, 50492, 62312,
                                                                       35099, 35162, 51052,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 63104, 0, 3,
                                                                       62312, 50520, 62348,
                                                                       35162, 35225, 51136,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 63212, 0, 3,
                                                                       62348, 50548, 62384,
                                                                       35225, 35288, 51220,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 63320, 0, 3,
                                                                       62420, 50604, 62456,
                                                                       35414, 35477, 51304,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 63428, 0, 3,
                                                                       62456, 50632, 62492,
                                                                       35477, 35540, 51388,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 63536, 0, 3,
                                                                       62492, 50660, 62528,
                                                                       35540, 35603, 51472,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 63644, 0, 3,
                                                                       62528, 50688, 62564,
                                                                       35603, 35666, 51556,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 63752, 0, 3,
                                                                       62564, 50716, 62600,
                                                                       35666, 35729, 51640,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 63860, 0, 3,
                                                                       62600, 50744, 62636,
                                                                       35729, 35792, 51724,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 63968, 0, 3,
                                                                       62672, 50800, 62780,
                                                                       35918, 36044, 51808,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 64184, 0, 3,
                                                                       62780, 50884, 62888,
                                                                       36044, 36170, 51976,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 64400, 0, 3,
                                                                       62888, 50968, 62996,
                                                                       36170, 36296, 52144,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 64616, 0, 3,
                                                                       62996, 51052, 63104,
                                                                       36296, 36422, 52312,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 64832, 0, 3,
                                                                       63104, 51136, 63212,
                                                                       36422, 36548, 52480,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 65048, 0, 3,
                                                                       63320, 51304, 63428,
                                                                       36800, 36926, 52648,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 65264, 0, 3,
                                                                       63428, 51388, 63536,
                                                                       36926, 37052, 52816,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 65480, 0, 3,
                                                                       63536, 51472, 63644,
                                                                       37052, 37178, 52984,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 65696, 0, 3,
                                                                       63644, 51556, 63752,
                                                                       37178, 37304, 53152,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 65912, 0, 3,
                                                                       63752, 51640, 63860,
                                                                       37304, 37430, 53320,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 66128, 0, 3,
                                                                       63968, 51808, 64184,
                                                                       37682, 37892, 53488,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 66488, 0, 3,
                                                                       64184, 51976, 64400,
                                                                       37892, 38102, 53768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 66848, 0, 3,
                                                                       64400, 52144, 64616,
                                                                       38102, 38312, 54048,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 67208, 0, 3,
                                                                       64616, 52312, 64832,
                                                                       38312, 38522, 54328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 67568, 0, 3,
                                                                       65048, 52648, 65264,
                                                                       38942, 39152, 54608,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 67928, 0, 3,
                                                                       65264, 52816, 65480,
                                                                       39152, 39362, 54888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 68288, 0, 3,
                                                                       65480, 52984, 65696,
                                                                       39362, 39572, 55168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 68648, 0, 3,
                                                                       65696, 53152, 65912,
                                                                       39572, 39782, 55448,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 69008, 0, 3,
                                                                       66128, 53488, 66488,
                                                                       40202, 40517, 55728,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 69548, 0, 3,
                                                                       66488, 53768, 66848,
                                                                       40517, 40832, 56148,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 70088, 0, 3,
                                                                       66848, 54048, 67208,
                                                                       40832, 41147, 56568,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 70628, 0, 3,
                                                                       67568, 54608, 67928,
                                                                       41777, 42092, 56988,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 71168, 0, 3,
                                                                       67928, 54888, 68288,
                                                                       42092, 42407, 57408,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 71708, 0, 3,
                                                                       68288, 55168, 68648,
                                                                       42407, 42722, 57828,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 72248, 0, 3,
                                                                       69008, 55728, 69548,
                                                                       43352, 43793, 58248,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 73004, 0, 3,
                                                                       69548, 56148, 70088,
                                                                       43793, 44234, 58836,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 73760, 0, 3,
                                                                       70628, 56988, 71168,
                                                                       45116, 45557, 59424,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 74516, 0, 3,
                                                                       71168, 57408, 71708,
                                                                       45557, 45998, 60012,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 75272, 0, 3,
                                                                       72248, 58248, 73004,
                                                                       46880, 47468, 60600,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 76280, 0, 3,
                                                                       73760, 59424, 74516,
                                                                       48644, 49232, 61384,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 77288, 76280, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 78296, 75272, 1008, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 79304, 77288, 28, 1, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 79304, 15, nmax);

        simdtrf::transform_k_inner(buffer, 79304, 78296, 28, 1, nmax);

        simdtrf::transform_i_outer(values + 195 * nvalues + n * npairs, nvalues, buffer, 79304,
                                   15, nmax);
    }

    for (size_t m = 0; m < 390; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
