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


#include "SimdThreeCenterElectronRepulsionRsRecIGD.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
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
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_igd_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_igd_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 49552, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1170 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 49552, 21172, 4730, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 12,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 20, 3, 12,
                                                             ncols, fj, i * nprim_b + j, fq);

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

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1436, 0, 3, 708,
                                                                       729, 1044, 1072, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1472, 0, 3, 729,
                                                                       750, 1072, 1100, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1508, 0, 3, 750,
                                                                       771, 1100, 1128, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1544, 0, 3, 771,
                                                                       792, 1128, 1156, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1580, 0, 3, 792,
                                                                       813, 1156, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1616, 0, 3, 813,
                                                                       834, 1184, 1212, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1652, 0, 3, 876,
                                                                       897, 1240, 1268, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 897,
                                                                       918, 1268, 1296, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1724, 0, 3, 918,
                                                                       939, 1296, 1324, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1760, 0, 3, 939,
                                                                       960, 1324, 1352, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1796, 0, 3, 960,
                                                                       981, 1352, 1380, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1832, 0, 3, 981,
                                                                       1002, 1380, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1868, 0, 3, 1044,
                                                                       1072, 1436, 1472, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1913, 0, 3, 1072,
                                                                       1100, 1472, 1508, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1958, 0, 3, 1100,
                                                                       1128, 1508, 1544, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2003, 0, 3, 1128,
                                                                       1156, 1544, 1580, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2048, 0, 3, 1156,
                                                                       1184, 1580, 1616, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2093, 0, 3, 1240,
                                                                       1268, 1652, 1688, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2138, 0, 3, 1268,
                                                                       1296, 1688, 1724, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2183, 0, 3, 1296,
                                                                       1324, 1724, 1760, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2228, 0, 3, 1324,
                                                                       1352, 1760, 1796, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2273, 0, 3, 1352,
                                                                       1380, 1796, 1832, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2318, 0, 3, 1436,
                                                                       1472, 1868, 1913, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2373, 0, 3, 1472,
                                                                       1508, 1913, 1958, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2428, 0, 3, 1508,
                                                                       1544, 1958, 2003, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2483, 0, 3, 1544,
                                                                       1580, 2003, 2048, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2538, 0, 3, 1652,
                                                                       1688, 2093, 2138, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2593, 0, 3, 1688,
                                                                       1724, 2138, 2183, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2648, 0, 3, 1724,
                                                                       1760, 2183, 2228, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2703, 0, 3, 1760,
                                                                       1796, 2228, 2273, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2758, 0, 3, 1868,
                                                                       1913, 2318, 2373, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2824, 0, 3, 1913,
                                                                       1958, 2373, 2428, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2890, 0, 3, 1958,
                                                                       2003, 2428, 2483, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2956, 0, 3, 2093,
                                                                       2138, 2538, 2593, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3022, 0, 3, 2138,
                                                                       2183, 2593, 2648, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3088, 0, 3, 2183,
                                                                       2228, 2648, 2703, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3154, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3157, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3160, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3163, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3166, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3169, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3172, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3175, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3178, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3181, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3184, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3187, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3190, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3193, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3196, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3199, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3202, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3205, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3208, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3211, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3214, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3217, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3220, 3, 9, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3229, 3, 10, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3238, 3, 11, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3247, 3, 12, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3256, 3, 13, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3265, 3, 14, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3274, 3, 15, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3283, 3, 16, 61,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3292, 3, 17, 64,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3301, 3, 18, 67,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3310, 3, 23, 76,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3319, 3, 24, 79,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3328, 3, 25, 82,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3337, 3, 26, 85,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3346, 3, 27, 88,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3355, 3, 28, 91,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3364, 3, 29, 94,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3373, 3, 30, 97,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3382, 3, 31, 100,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3391, 3, 32, 103,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3400, 3, 40, 118,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3418, 3, 43, 124,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3436, 3, 46, 130,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3454, 3, 49, 136,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3472, 3, 52, 142,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3490, 3, 55, 148,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3508, 3, 58, 154,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3526, 3, 61, 160,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3544, 3, 64, 166,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3562, 3, 76, 184,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3580, 3, 79, 190,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3598, 3, 82, 196,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3616, 3, 85, 202,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3634, 3, 88, 208,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3652, 3, 91, 214,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3670, 3, 94, 220,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3688, 3, 97, 226,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3706, 3, 100, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3724, 3, 118, 258,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3754, 3, 124, 268,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3784, 3, 130, 278,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3814, 3, 136, 288,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3844, 3, 142, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3874, 3, 148, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3904, 3, 154, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3934, 3, 160, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3964, 3, 184, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3994, 3, 190, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4024, 3, 196, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4054, 3, 202, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4084, 3, 208, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4114, 3, 214, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4144, 3, 220, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4174, 3, 226, 428,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4204, 3, 258, 468,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4249, 3, 268, 483,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4294, 3, 278, 498,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4339, 3, 288, 513,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4384, 3, 298, 528,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4429, 3, 308, 543,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4474, 3, 318, 558,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4519, 3, 358, 603,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4564, 3, 368, 618,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4609, 3, 378, 633,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4654, 3, 388, 648,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4699, 3, 398, 663,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4744, 3, 408, 678,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4789, 3, 418, 693,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4834, 3, 468, 750,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4897, 3, 483, 771,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4960, 3, 498, 792,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5023, 3, 513, 813,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5086, 3, 528, 834,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5149, 3, 543, 855,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5212, 3, 603, 918,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5275, 3, 618, 939,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5338, 3, 633, 960,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5401, 3, 648, 981,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5464, 3, 663,
                                                                       1002, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5527, 3, 678,
                                                                       1023, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5590, 3, 750,
                                                                       1100, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5674, 3, 771,
                                                                       1128, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5758, 3, 792,
                                                                       1156, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5842, 3, 813,
                                                                       1184, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5926, 3, 834,
                                                                       1212, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6010, 3, 918,
                                                                       1296, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6094, 3, 939,
                                                                       1324, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6178, 3, 960,
                                                                       1352, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6262, 3, 981,
                                                                       1380, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6346, 3, 1002,
                                                                       1408, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6430, 3, 1100,
                                                                       1508, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6538, 3, 1128,
                                                                       1544, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6646, 3, 1156,
                                                                       1580, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6754, 3, 1184,
                                                                       1616, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6862, 3, 1296,
                                                                       1724, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6970, 3, 1324,
                                                                       1760, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7078, 3, 1352,
                                                                       1796, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7186, 3, 1380,
                                                                       1832, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7294, 3, 1508,
                                                                       1958, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7429, 3, 1544,
                                                                       2003, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7564, 3, 1580,
                                                                       2048, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7699, 3, 1724,
                                                                       2183, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7834, 3, 1760,
                                                                       2228, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7969, 3, 1796,
                                                                       2273, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8104, 3, 1958,
                                                                       2428, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8269, 3, 2003,
                                                                       2483, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8434, 3, 2183,
                                                                       2648, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8599, 3, 2228,
                                                                       2703, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 8764, 3, 2428,
                                                                       2890, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 8962, 3, 2648,
                                                                       3088, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9160, 3, 7, 8,
                                                                       3154, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9166, 3, 8, 9,
                                                                       3157, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9172, 3, 9, 10,
                                                                       3160, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9178, 3, 10, 11,
                                                                       3163, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9184, 3, 11, 12,
                                                                       3166, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9190, 3, 12, 13,
                                                                       3169, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9196, 3, 13, 14,
                                                                       3172, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9202, 3, 14, 15,
                                                                       3175, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9208, 3, 15, 16,
                                                                       3178, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9214, 3, 16, 17,
                                                                       3181, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9220, 3, 17, 18,
                                                                       3184, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9226, 3, 21, 22,
                                                                       3187, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9232, 3, 22, 23,
                                                                       3190, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9238, 3, 23, 24,
                                                                       3193, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9244, 3, 24, 25,
                                                                       3196, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9250, 3, 25, 26,
                                                                       3199, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9256, 3, 26, 27,
                                                                       3202, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9262, 3, 27, 28,
                                                                       3205, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9268, 3, 28, 29,
                                                                       3208, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9274, 3, 29, 30,
                                                                       3211, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9280, 3, 30, 31,
                                                                       3214, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9286, 3, 31, 32,
                                                                       3217, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9292, 0, 3, 9160,
                                                                       3154, 9166, 3220, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9310, 0, 3, 9166,
                                                                       3157, 9172, 3229, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9328, 0, 3, 9172,
                                                                       3160, 9178, 3238, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9346, 0, 3, 9178,
                                                                       3163, 9184, 3247, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9364, 0, 3, 9184,
                                                                       3166, 9190, 3256, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9382, 0, 3, 9190,
                                                                       3169, 9196, 3265, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9400, 0, 3, 9196,
                                                                       3172, 9202, 3274, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9418, 0, 3, 9202,
                                                                       3175, 9208, 3283, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9436, 0, 3, 9208,
                                                                       3178, 9214, 3292, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9454, 0, 3, 9214,
                                                                       3181, 9220, 3301, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9472, 0, 3, 9226,
                                                                       3187, 9232, 3310, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9490, 0, 3, 9232,
                                                                       3190, 9238, 3319, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9508, 0, 3, 9238,
                                                                       3193, 9244, 3328, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9526, 0, 3, 9244,
                                                                       3196, 9250, 3337, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9544, 0, 3, 9250,
                                                                       3199, 9256, 3346, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9562, 0, 3, 9256,
                                                                       3202, 9262, 3355, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9580, 0, 3, 9262,
                                                                       3205, 9268, 3364, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9598, 0, 3, 9268,
                                                                       3208, 9274, 3373, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9616, 0, 3, 9274,
                                                                       3211, 9280, 3382, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9634, 0, 3, 9280,
                                                                       3214, 9286, 3391, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9652, 0, 3, 9292,
                                                                       3220, 9310, 106, 112,
                                                                       3400, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9688, 0, 3, 9310,
                                                                       3229, 9328, 112, 118,
                                                                       3418, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9724, 0, 3, 9328,
                                                                       3238, 9346, 118, 124,
                                                                       3436, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9760, 0, 3, 9346,
                                                                       3247, 9364, 124, 130,
                                                                       3454, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9796, 0, 3, 9364,
                                                                       3256, 9382, 130, 136,
                                                                       3472, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9832, 0, 3, 9382,
                                                                       3265, 9400, 136, 142,
                                                                       3490, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9868, 0, 3, 9400,
                                                                       3274, 9418, 142, 148,
                                                                       3508, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9904, 0, 3, 9418,
                                                                       3283, 9436, 148, 154,
                                                                       3526, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9940, 0, 3, 9436,
                                                                       3292, 9454, 154, 160,
                                                                       3544, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9976, 0, 3, 9472,
                                                                       3310, 9490, 172, 178,
                                                                       3562, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10012, 0, 3, 9490,
                                                                       3319, 9508, 178, 184,
                                                                       3580, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10048, 0, 3, 9508,
                                                                       3328, 9526, 184, 190,
                                                                       3598, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10084, 0, 3, 9526,
                                                                       3337, 9544, 190, 196,
                                                                       3616, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10120, 0, 3, 9544,
                                                                       3346, 9562, 196, 202,
                                                                       3634, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10156, 0, 3, 9562,
                                                                       3355, 9580, 202, 208,
                                                                       3652, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10192, 0, 3, 9580,
                                                                       3364, 9598, 208, 214,
                                                                       3670, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10228, 0, 3, 9598,
                                                                       3373, 9616, 214, 220,
                                                                       3688, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10264, 0, 3, 9616,
                                                                       3382, 9634, 220, 226,
                                                                       3706, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10300, 0, 3, 9652,
                                                                       3400, 9688, 238, 248,
                                                                       3724, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10360, 0, 3, 9688,
                                                                       3418, 9724, 248, 258,
                                                                       3754, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10420, 0, 3, 9724,
                                                                       3436, 9760, 258, 268,
                                                                       3784, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10480, 0, 3, 9760,
                                                                       3454, 9796, 268, 278,
                                                                       3814, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10540, 0, 3, 9796,
                                                                       3472, 9832, 278, 288,
                                                                       3844, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10600, 0, 3, 9832,
                                                                       3490, 9868, 288, 298,
                                                                       3874, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10660, 0, 3, 9868,
                                                                       3508, 9904, 298, 308,
                                                                       3904, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10720, 0, 3, 9904,
                                                                       3526, 9940, 308, 318,
                                                                       3934, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10780, 0, 3, 9976,
                                                                       3562, 10012, 338, 348,
                                                                       3964, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10840, 0, 3,
                                                                       10012, 3580, 10048, 348,
                                                                       358, 3994, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10900, 0, 3,
                                                                       10048, 3598, 10084, 358,
                                                                       368, 4024, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10960, 0, 3,
                                                                       10084, 3616, 10120, 368,
                                                                       378, 4054, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11020, 0, 3,
                                                                       10120, 3634, 10156, 378,
                                                                       388, 4084, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11080, 0, 3,
                                                                       10156, 3652, 10192, 388,
                                                                       398, 4114, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11140, 0, 3,
                                                                       10192, 3670, 10228, 398,
                                                                       408, 4144, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11200, 0, 3,
                                                                       10228, 3688, 10264, 408,
                                                                       418, 4174, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11260, 0, 3,
                                                                       10300, 3724, 10360, 438,
                                                                       453, 4204, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11350, 0, 3,
                                                                       10360, 3754, 10420, 453,
                                                                       468, 4249, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11440, 0, 3,
                                                                       10420, 3784, 10480, 468,
                                                                       483, 4294, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11530, 0, 3,
                                                                       10480, 3814, 10540, 483,
                                                                       498, 4339, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11620, 0, 3,
                                                                       10540, 3844, 10600, 498,
                                                                       513, 4384, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11710, 0, 3,
                                                                       10600, 3874, 10660, 513,
                                                                       528, 4429, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11800, 0, 3,
                                                                       10660, 3904, 10720, 528,
                                                                       543, 4474, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11890, 0, 3,
                                                                       10780, 3964, 10840, 573,
                                                                       588, 4519, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11980, 0, 3,
                                                                       10840, 3994, 10900, 588,
                                                                       603, 4564, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12070, 0, 3,
                                                                       10900, 4024, 10960, 603,
                                                                       618, 4609, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12160, 0, 3,
                                                                       10960, 4054, 11020, 618,
                                                                       633, 4654, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12250, 0, 3,
                                                                       11020, 4084, 11080, 633,
                                                                       648, 4699, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12340, 0, 3,
                                                                       11080, 4114, 11140, 648,
                                                                       663, 4744, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12430, 0, 3,
                                                                       11140, 4144, 11200, 663,
                                                                       678, 4789, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12520, 0, 3,
                                                                       11260, 4204, 11350, 708,
                                                                       729, 4834, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12646, 0, 3,
                                                                       11350, 4249, 11440, 729,
                                                                       750, 4897, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12772, 0, 3,
                                                                       11440, 4294, 11530, 750,
                                                                       771, 4960, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12898, 0, 3,
                                                                       11530, 4339, 11620, 771,
                                                                       792, 5023, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13024, 0, 3,
                                                                       11620, 4384, 11710, 792,
                                                                       813, 5086, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13150, 0, 3,
                                                                       11710, 4429, 11800, 813,
                                                                       834, 5149, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13276, 0, 3,
                                                                       11890, 4519, 11980, 876,
                                                                       897, 5212, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13402, 0, 3,
                                                                       11980, 4564, 12070, 897,
                                                                       918, 5275, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13528, 0, 3,
                                                                       12070, 4609, 12160, 918,
                                                                       939, 5338, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13654, 0, 3,
                                                                       12160, 4654, 12250, 939,
                                                                       960, 5401, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13780, 0, 3,
                                                                       12250, 4699, 12340, 960,
                                                                       981, 5464, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13906, 0, 3,
                                                                       12340, 4744, 12430, 981,
                                                                       1002, 5527, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14032, 0, 3,
                                                                       12520, 4834, 12646, 1044,
                                                                       1072, 5590, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14200, 0, 3,
                                                                       12646, 4897, 12772, 1072,
                                                                       1100, 5674, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14368, 0, 3,
                                                                       12772, 4960, 12898, 1100,
                                                                       1128, 5758, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14536, 0, 3,
                                                                       12898, 5023, 13024, 1128,
                                                                       1156, 5842, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14704, 0, 3,
                                                                       13024, 5086, 13150, 1156,
                                                                       1184, 5926, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14872, 0, 3,
                                                                       13276, 5212, 13402, 1240,
                                                                       1268, 6010, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15040, 0, 3,
                                                                       13402, 5275, 13528, 1268,
                                                                       1296, 6094, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15208, 0, 3,
                                                                       13528, 5338, 13654, 1296,
                                                                       1324, 6178, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15376, 0, 3,
                                                                       13654, 5401, 13780, 1324,
                                                                       1352, 6262, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15544, 0, 3,
                                                                       13780, 5464, 13906, 1352,
                                                                       1380, 6346, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15712, 0, 3,
                                                                       14032, 5590, 14200, 1436,
                                                                       1472, 6430, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15928, 0, 3,
                                                                       14200, 5674, 14368, 1472,
                                                                       1508, 6538, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16144, 0, 3,
                                                                       14368, 5758, 14536, 1508,
                                                                       1544, 6646, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16360, 0, 3,
                                                                       14536, 5842, 14704, 1544,
                                                                       1580, 6754, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16576, 0, 3,
                                                                       14872, 6010, 15040, 1652,
                                                                       1688, 6862, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16792, 0, 3,
                                                                       15040, 6094, 15208, 1688,
                                                                       1724, 6970, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17008, 0, 3,
                                                                       15208, 6178, 15376, 1724,
                                                                       1760, 7078, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17224, 0, 3,
                                                                       15376, 6262, 15544, 1760,
                                                                       1796, 7186, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17440, 0, 3,
                                                                       15712, 6430, 15928, 1868,
                                                                       1913, 7294, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17710, 0, 3,
                                                                       15928, 6538, 16144, 1913,
                                                                       1958, 7429, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17980, 0, 3,
                                                                       16144, 6646, 16360, 1958,
                                                                       2003, 7564, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18250, 0, 3,
                                                                       16576, 6862, 16792, 2093,
                                                                       2138, 7699, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18520, 0, 3,
                                                                       16792, 6970, 17008, 2138,
                                                                       2183, 7834, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18790, 0, 3,
                                                                       17008, 7078, 17224, 2183,
                                                                       2228, 7969, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 19060, 0, 3,
                                                                       17440, 7294, 17710, 2318,
                                                                       2373, 8104, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 19390, 0, 3,
                                                                       17710, 7429, 17980, 2373,
                                                                       2428, 8269, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 19720, 0, 3,
                                                                       18250, 7699, 18520, 2538,
                                                                       2593, 8434, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 20050, 0, 3,
                                                                       18520, 7834, 18790, 2593,
                                                                       2648, 8599, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 20380, 0, 3,
                                                                       19060, 8104, 19390, 2758,
                                                                       2824, 8764, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 20776, 0, 3,
                                                                       19720, 8434, 20050, 2956,
                                                                       3022, 8962, ncols, gamma,
                                                                       p, q);

                    simdfunc::contract_primitives(buffer, 21172, 14032, 168, ncols);

                    simdfunc::contract_primitives(buffer, 21480, 14872, 168, ncols);

                    simdfunc::contract_primitives(buffer, 21788, 15712, 216, ncols);

                    simdfunc::contract_primitives(buffer, 22184, 16576, 216, ncols);

                    simdfunc::contract_primitives(buffer, 22580, 17440, 270, ncols);

                    simdfunc::contract_primitives(buffer, 23075, 18250, 270, ncols);

                    simdfunc::contract_primitives(buffer, 23570, 19060, 330, ncols);

                    simdfunc::contract_primitives(buffer, 24175, 19720, 330, ncols);

                    simdfunc::contract_primitives(buffer, 24780, 20380, 396, ncols);

                    simdfunc::contract_primitives(buffer, 25506, 20776, 396, ncols);
                }
            }
        }

        simdtrf::transform_d_inner(buffer, 21340, 21172, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 21648, 21480, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 22004, 21788, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 22400, 22184, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 22850, 22580, 45, 1, nmax);

        simdtrf::transform_d_inner(buffer, 23345, 23075, 45, 1, nmax);

        simdtrf::transform_d_inner(buffer, 23900, 23570, 55, 1, nmax);

        simdtrf::transform_d_inner(buffer, 24505, 24175, 55, 1, nmax);

        simdtrf::transform_d_inner(buffer, 25176, 24780, 66, 1, nmax);

        simdtrf::transform_d_inner(buffer, 25902, 25506, 66, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 26232, 21340, 22004, 5, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 26652, 21648, 22400, 5, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 27072, 22004, 22850, 5, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 27612, 22400, 23345, 5, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 28152, 22850, 23900, 5, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 28827, 23345, 24505, 5, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 29502, 23900, 25176, 5, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 30327, 24505, 25902, 5, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 31152, 26232, 27072, 5, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 31992, 26652, 27612, 5, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 32832, 27072, 28152, 5, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 33912, 27612, 28827, 5, nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 34992, 28152, 29502, 5, nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 36342, 28827, 30327, 5, nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 37692, 31152, 32832, 5, nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 39092, 31992, 33912, 5, nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 40492, 32832, 34992, 5, nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 42292, 33912, 36342, 5, nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 44092, 37692, 40492, 5, nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 46192, 39092, 42292, 5, nmax);

        simdtrf::transform_g_inner(buffer, 48292, 46192, 28, 5, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 48292, 45, nmax);

        simdtrf::transform_g_inner(buffer, 48292, 44092, 28, 5, nmax);

        simdtrf::transform_i_outer(values + 585 * nvalues + n * npairs, nvalues, buffer, 48292,
                                   45, nmax);
    }

    for (size_t m = 0; m < 1170; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
