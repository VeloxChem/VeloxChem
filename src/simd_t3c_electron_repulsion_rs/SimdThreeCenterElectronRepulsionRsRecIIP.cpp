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


#include "SimdThreeCenterElectronRepulsionRsRecIIP.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSS.hpp"
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
#include "SimdTransformI.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_iip_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_iip_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 58164, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1014 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 58164, 6042, 4515, dimensions);

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

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3154, 0, 3, 2318,
                                                                       2373, 2758, 2824, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3232, 0, 3, 2373,
                                                                       2428, 2824, 2890, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3310, 0, 3, 2538,
                                                                       2593, 2956, 3022, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3388, 0, 3, 2593,
                                                                       2648, 3022, 3088, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 3466, 0, 3, 2758,
                                                                       2824, 3154, 3232, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 3557, 0, 3, 2956,
                                                                       3022, 3310, 3388, ncols,
                                                                       gamma, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3648, 3, 708,
                                                                       1044, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3732, 3, 876,
                                                                       1240, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3816, 3, 1044,
                                                                       1436, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3924, 3, 1240,
                                                                       1652, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 4032, 3, 1436,
                                                                       1868, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 4167, 3, 1652,
                                                                       2093, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 4302, 3, 1868,
                                                                       2318, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 4467, 3, 2093,
                                                                       2538, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 4632, 3, 2318,
                                                                       2758, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 4830, 3, 2538,
                                                                       2956, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 5028, 3, 2758,
                                                                       3154, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 5262, 3, 2956,
                                                                       3310, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 5496, 3, 3154,
                                                                       3466, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 5769, 3, 3310,
                                                                       3557, ncols, p, q);

                    simdfunc::contract_primitives(buffer, 6042, 3648, 84, ncols);

                    simdfunc::contract_primitives(buffer, 6210, 3732, 84, ncols);

                    simdfunc::contract_primitives(buffer, 6378, 3816, 108, ncols);

                    simdfunc::contract_primitives(buffer, 6594, 3924, 108, ncols);

                    simdfunc::contract_primitives(buffer, 6810, 4032, 135, ncols);

                    simdfunc::contract_primitives(buffer, 7080, 4167, 135, ncols);

                    simdfunc::contract_primitives(buffer, 7350, 4302, 165, ncols);

                    simdfunc::contract_primitives(buffer, 7680, 4467, 165, ncols);

                    simdfunc::contract_primitives(buffer, 8010, 4632, 198, ncols);

                    simdfunc::contract_primitives(buffer, 8406, 4830, 198, ncols);

                    simdfunc::contract_primitives(buffer, 8802, 5028, 234, ncols);

                    simdfunc::contract_primitives(buffer, 9270, 5262, 234, ncols);

                    simdfunc::contract_primitives(buffer, 9738, 5496, 273, ncols);

                    simdfunc::contract_primitives(buffer, 10284, 5769, 273, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 6126, 6042, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6294, 6210, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6486, 6378, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6702, 6594, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6945, 6810, 45, 1, nmax);

        simdtrf::transform_p_inner(buffer, 7215, 7080, 45, 1, nmax);

        simdtrf::transform_p_inner(buffer, 7515, 7350, 55, 1, nmax);

        simdtrf::transform_p_inner(buffer, 7845, 7680, 55, 1, nmax);

        simdtrf::transform_p_inner(buffer, 8208, 8010, 66, 1, nmax);

        simdtrf::transform_p_inner(buffer, 8604, 8406, 66, 1, nmax);

        simdtrf::transform_p_inner(buffer, 9036, 8802, 78, 1, nmax);

        simdtrf::transform_p_inner(buffer, 9504, 9270, 78, 1, nmax);

        simdtrf::transform_p_inner(buffer, 10011, 9738, 91, 1, nmax);

        simdtrf::transform_p_inner(buffer, 10557, 10284, 91, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 10830, 6126, 6486, 3, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 11082, 6294, 6702, 3, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 11334, 6486, 6945, 3, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 11658, 6702, 7215, 3, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 11982, 6945, 7515, 3, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 12387, 7215, 7845, 3, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 12792, 7515, 8208, 3, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 13287, 7845, 8604, 3, nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 13782, 8208, 9036, 3, nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 14376, 8604, 9504, 3, nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 14970, 9036, 10011, 3, nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 15672, 9504, 10557, 3, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 16374, 10830, 11334, 3, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 16878, 11082, 11658, 3, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 17382, 11334, 11982, 3, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 18030, 11658, 12387, 3, nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 18678, 11982, 12792, 3, nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 19488, 12387, 13287, 3, nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 20298, 12792, 13782, 3, nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 21288, 13287, 14376, 3, nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 22278, 13782, 14970, 3, nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 23466, 14376, 15672, 3, nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 24654, 16374, 17382, 3, nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 25494, 16878, 18030, 3, nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 26334, 17382, 18678, 3, nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 27414, 18030, 19488, 3, nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 28494, 18678, 20298, 3, nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 29844, 19488, 21288, 3, nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 31194, 20298, 22278, 3, nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 32844, 21288, 23466, 3, nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 34494, 24654, 26334, 3, nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 35754, 25494, 27414, 3, nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 37014, 26334, 28494, 3, nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 38634, 27414, 29844, 3, nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 40254, 28494, 31194, 3, nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 42279, 29844, 32844, 3, nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 44304, 34494, 37014, 3, nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 46068, 35754, 38634, 3, nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 47832, 37014, 40254, 3, nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 50100, 38634, 42279, 3, nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 52368, 44304, 47832, 3, nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 54720, 46068, 50100, 3, nmax);

        simdtrf::transform_i_inner(buffer, 57072, 54720, 28, 3, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 57072, 39, nmax);

        simdtrf::transform_i_inner(buffer, 57072, 52368, 28, 3, nmax);

        simdtrf::transform_i_outer(values + 507 * nvalues + n * npairs, nvalues, buffer, 57072,
                                   39, nmax);
    }

    for (size_t m = 0; m < 1014; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
