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


#include "SimdThreeCenterElectronRepulsionRsRecIFF.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIF.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_iff_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_iff_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 52950, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1274 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 52950, 32128, 5191, dimensions);

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
                                                            4, 5, 6, 7, 8, 9, 10, 11, 12}, ncols,
                                                            fj, i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 19, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12}, ncols, fj,
                                                        i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 32, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 35, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 38, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 41, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 44, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 47, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 50, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 53, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 56, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 59, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 62, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 65, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 68, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 71, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 74, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 77, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 80, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 83, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 86, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 89, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 92, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 95, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 98, 0, 3, 7, 8,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 104, 0, 3, 8, 9,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 110, 0, 3, 9, 10,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 116, 0, 3, 10, 11,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 122, 0, 3, 11, 12,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 128, 0, 3, 12, 13,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 134, 0, 3, 13, 14,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 140, 0, 3, 14, 15,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 146, 0, 3, 15, 16,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 152, 0, 3, 16, 17,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 158, 0, 3, 20, 21,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 164, 0, 3, 21, 22,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 170, 0, 3, 22, 23,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 176, 0, 3, 23, 24,
                                                                       74, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 182, 0, 3, 24, 25,
                                                                       77, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 188, 0, 3, 25, 26,
                                                                       80, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 194, 0, 3, 26, 27,
                                                                       83, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 200, 0, 3, 27, 28,
                                                                       86, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 206, 0, 3, 28, 29,
                                                                       89, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 212, 0, 3, 29, 30,
                                                                       92, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 218, 0, 3, 32, 35,
                                                                       98, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 228, 0, 3, 35, 38,
                                                                       104, 110, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 238, 0, 3, 38, 41,
                                                                       110, 116, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 248, 0, 3, 41, 44,
                                                                       116, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 258, 0, 3, 44, 47,
                                                                       122, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 268, 0, 3, 47, 50,
                                                                       128, 134, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 278, 0, 3, 50, 53,
                                                                       134, 140, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 288, 0, 3, 53, 56,
                                                                       140, 146, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 298, 0, 3, 56, 59,
                                                                       146, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 308, 0, 3, 65, 68,
                                                                       158, 164, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 318, 0, 3, 68, 71,
                                                                       164, 170, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 328, 0, 3, 71, 74,
                                                                       170, 176, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 338, 0, 3, 74, 77,
                                                                       176, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 348, 0, 3, 77, 80,
                                                                       182, 188, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 80, 83,
                                                                       188, 194, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 368, 0, 3, 83, 86,
                                                                       194, 200, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 378, 0, 3, 86, 89,
                                                                       200, 206, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 388, 0, 3, 89, 92,
                                                                       206, 212, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 398, 0, 3, 98,
                                                                       104, 218, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 413, 0, 3, 104,
                                                                       110, 228, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 110,
                                                                       116, 238, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 443, 0, 3, 116,
                                                                       122, 248, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 458, 0, 3, 122,
                                                                       128, 258, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 473, 0, 3, 128,
                                                                       134, 268, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 488, 0, 3, 134,
                                                                       140, 278, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 503, 0, 3, 140,
                                                                       146, 288, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 518, 0, 3, 158,
                                                                       164, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 533, 0, 3, 164,
                                                                       170, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 548, 0, 3, 170,
                                                                       176, 328, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 563, 0, 3, 176,
                                                                       182, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 578, 0, 3, 182,
                                                                       188, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 593, 0, 3, 188,
                                                                       194, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 608, 0, 3, 194,
                                                                       200, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 623, 0, 3, 200,
                                                                       206, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 638, 0, 3, 218,
                                                                       228, 398, 413, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 659, 0, 3, 228,
                                                                       238, 413, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 680, 0, 3, 238,
                                                                       248, 428, 443, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 701, 0, 3, 248,
                                                                       258, 443, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 722, 0, 3, 258,
                                                                       268, 458, 473, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 743, 0, 3, 268,
                                                                       278, 473, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 764, 0, 3, 278,
                                                                       288, 488, 503, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 785, 0, 3, 308,
                                                                       318, 518, 533, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 806, 0, 3, 318,
                                                                       328, 533, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 827, 0, 3, 328,
                                                                       338, 548, 563, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 848, 0, 3, 338,
                                                                       348, 563, 578, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 869, 0, 3, 348,
                                                                       358, 578, 593, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 890, 0, 3, 358,
                                                                       368, 593, 608, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 911, 0, 3, 368,
                                                                       378, 608, 623, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 932, 0, 3, 398,
                                                                       413, 638, 659, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 960, 0, 3, 413,
                                                                       428, 659, 680, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 988, 0, 3, 428,
                                                                       443, 680, 701, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1016, 0, 3, 443,
                                                                       458, 701, 722, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1044, 0, 3, 458,
                                                                       473, 722, 743, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1072, 0, 3, 473,
                                                                       488, 743, 764, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1100, 0, 3, 518,
                                                                       533, 785, 806, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1128, 0, 3, 533,
                                                                       548, 806, 827, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1156, 0, 3, 548,
                                                                       563, 827, 848, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1184, 0, 3, 563,
                                                                       578, 848, 869, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1212, 0, 3, 578,
                                                                       593, 869, 890, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1240, 0, 3, 593,
                                                                       608, 890, 911, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 638,
                                                                       659, 932, 960, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1304, 0, 3, 659,
                                                                       680, 960, 988, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1340, 0, 3, 680,
                                                                       701, 988, 1016, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1376, 0, 3, 701,
                                                                       722, 1016, 1044, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1412, 0, 3, 722,
                                                                       743, 1044, 1072, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1448, 0, 3, 785,
                                                                       806, 1100, 1128, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1484, 0, 3, 806,
                                                                       827, 1128, 1156, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 827,
                                                                       848, 1156, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1556, 0, 3, 848,
                                                                       869, 1184, 1212, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1592, 0, 3, 869,
                                                                       890, 1212, 1240, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1628, 0, 3, 932,
                                                                       960, 1268, 1304, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1673, 0, 3, 960,
                                                                       988, 1304, 1340, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1718, 0, 3, 988,
                                                                       1016, 1340, 1376, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1763, 0, 3, 1016,
                                                                       1044, 1376, 1412, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1808, 0, 3, 1100,
                                                                       1128, 1448, 1484, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1853, 0, 3, 1128,
                                                                       1156, 1484, 1520, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1898, 0, 3, 1156,
                                                                       1184, 1520, 1556, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1943, 0, 3, 1184,
                                                                       1212, 1556, 1592, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1988, 0, 3, 1268,
                                                                       1304, 1628, 1673, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2043, 0, 3, 1304,
                                                                       1340, 1673, 1718, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2098, 0, 3, 1340,
                                                                       1376, 1718, 1763, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2153, 0, 3, 1448,
                                                                       1484, 1808, 1853, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2208, 0, 3, 1484,
                                                                       1520, 1853, 1898, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2263, 0, 3, 1520,
                                                                       1556, 1898, 1943, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2318, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2321, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2324, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2327, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2330, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2333, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2336, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2339, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2342, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2345, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2348, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2351, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2354, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2357, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2360, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2363, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2366, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2369, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2372, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2375, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2378, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2381, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2384, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2387, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2390, 3, 7, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2399, 3, 8, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2408, 3, 9, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2417, 3, 10, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2426, 3, 11, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2435, 3, 12, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2444, 3, 13, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2453, 3, 14, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2462, 3, 15, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2471, 3, 16, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2480, 3, 17, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2489, 3, 20, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2498, 3, 21, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2507, 3, 22, 71,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2516, 3, 23, 74,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2525, 3, 24, 77,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2534, 3, 25, 80,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2543, 3, 26, 83,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2552, 3, 27, 86,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2561, 3, 28, 89,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2570, 3, 29, 92,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2579, 3, 30, 95,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2588, 3, 32, 98,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2606, 3, 35, 104,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2624, 3, 38, 110,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2642, 3, 41, 116,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2660, 3, 44, 122,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2678, 3, 47, 128,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2696, 3, 50, 134,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2714, 3, 53, 140,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2732, 3, 56, 146,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2750, 3, 59, 152,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2768, 3, 65, 158,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2786, 3, 68, 164,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2804, 3, 71, 170,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2822, 3, 74, 176,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2840, 3, 77, 182,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2858, 3, 80, 188,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2876, 3, 83, 194,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2894, 3, 86, 200,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2912, 3, 89, 206,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2930, 3, 92, 212,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2948, 3, 98, 218,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2978, 3, 104, 228,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3008, 3, 110, 238,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3038, 3, 116, 248,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3068, 3, 122, 258,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3098, 3, 128, 268,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3128, 3, 134, 278,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3158, 3, 140, 288,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3188, 3, 146, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3218, 3, 158, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3248, 3, 164, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3278, 3, 170, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3308, 3, 176, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3338, 3, 182, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3368, 3, 188, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3398, 3, 194, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3428, 3, 200, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3458, 3, 206, 388,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3488, 3, 218, 398,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3533, 3, 228, 413,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3578, 3, 238, 428,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3623, 3, 248, 443,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3668, 3, 258, 458,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3713, 3, 268, 473,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3758, 3, 278, 488,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3803, 3, 288, 503,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3848, 3, 308, 518,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3893, 3, 318, 533,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3938, 3, 328, 548,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3983, 3, 338, 563,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4028, 3, 348, 578,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4073, 3, 358, 593,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4118, 3, 368, 608,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4163, 3, 378, 623,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4208, 3, 398, 638,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4271, 3, 413, 659,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4334, 3, 428, 680,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4397, 3, 443, 701,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4460, 3, 458, 722,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4523, 3, 473, 743,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4586, 3, 488, 764,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4649, 3, 518, 785,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4712, 3, 533, 806,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4775, 3, 548, 827,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4838, 3, 563, 848,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4901, 3, 578, 869,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4964, 3, 593, 890,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5027, 3, 608, 911,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5090, 3, 638, 932,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5174, 3, 659, 960,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5258, 3, 680, 988,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5342, 3, 701,
                                                                       1016, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5426, 3, 722,
                                                                       1044, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5510, 3, 743,
                                                                       1072, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5594, 3, 785,
                                                                       1100, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5678, 3, 806,
                                                                       1128, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5762, 3, 827,
                                                                       1156, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5846, 3, 848,
                                                                       1184, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5930, 3, 869,
                                                                       1212, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6014, 3, 890,
                                                                       1240, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6098, 3, 932,
                                                                       1268, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6206, 3, 960,
                                                                       1304, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6314, 3, 988,
                                                                       1340, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6422, 3, 1016,
                                                                       1376, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6530, 3, 1044,
                                                                       1412, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6638, 3, 1100,
                                                                       1448, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6746, 3, 1128,
                                                                       1484, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6854, 3, 1156,
                                                                       1520, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6962, 3, 1184,
                                                                       1556, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7070, 3, 1212,
                                                                       1592, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7178, 3, 1268,
                                                                       1628, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7313, 3, 1304,
                                                                       1673, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7448, 3, 1340,
                                                                       1718, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7583, 3, 1376,
                                                                       1763, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7718, 3, 1448,
                                                                       1808, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7853, 3, 1484,
                                                                       1853, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7988, 3, 1520,
                                                                       1898, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8123, 3, 1556,
                                                                       1943, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8258, 3, 1628,
                                                                       1988, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8423, 3, 1673,
                                                                       2043, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8588, 3, 1718,
                                                                       2098, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8753, 3, 1808,
                                                                       2153, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8918, 3, 1853,
                                                                       2208, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9083, 3, 1898,
                                                                       2263, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9248, 3, 7, 8,
                                                                       2324, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9254, 3, 8, 9,
                                                                       2327, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9260, 3, 9, 10,
                                                                       2330, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9266, 3, 10, 11,
                                                                       2333, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9272, 3, 11, 12,
                                                                       2336, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9278, 3, 12, 13,
                                                                       2339, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9284, 3, 13, 14,
                                                                       2342, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9290, 3, 14, 15,
                                                                       2345, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9296, 3, 15, 16,
                                                                       2348, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9302, 3, 16, 17,
                                                                       2351, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9308, 3, 20, 21,
                                                                       2360, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9314, 3, 21, 22,
                                                                       2363, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9320, 3, 22, 23,
                                                                       2366, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9326, 3, 23, 24,
                                                                       2369, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9332, 3, 24, 25,
                                                                       2372, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9338, 3, 25, 26,
                                                                       2375, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9344, 3, 26, 27,
                                                                       2378, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9350, 3, 27, 28,
                                                                       2381, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9356, 3, 28, 29,
                                                                       2384, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9362, 3, 29, 30,
                                                                       2387, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9368, 0, 3, 9248,
                                                                       2324, 9254, 2408, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9386, 0, 3, 9254,
                                                                       2327, 9260, 2417, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9404, 0, 3, 9260,
                                                                       2330, 9266, 2426, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9422, 0, 3, 9266,
                                                                       2333, 9272, 2435, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9440, 0, 3, 9272,
                                                                       2336, 9278, 2444, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9458, 0, 3, 9278,
                                                                       2339, 9284, 2453, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9476, 0, 3, 9284,
                                                                       2342, 9290, 2462, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9494, 0, 3, 9290,
                                                                       2345, 9296, 2471, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9512, 0, 3, 9296,
                                                                       2348, 9302, 2480, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9530, 0, 3, 9308,
                                                                       2360, 9314, 2507, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9548, 0, 3, 9314,
                                                                       2363, 9320, 2516, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9566, 0, 3, 9320,
                                                                       2366, 9326, 2525, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9584, 0, 3, 9326,
                                                                       2369, 9332, 2534, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9602, 0, 3, 9332,
                                                                       2372, 9338, 2543, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9620, 0, 3, 9338,
                                                                       2375, 9344, 2552, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9638, 0, 3, 9344,
                                                                       2378, 9350, 2561, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9656, 0, 3, 9350,
                                                                       2381, 9356, 2570, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9674, 0, 3, 9356,
                                                                       2384, 9362, 2579, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9692, 0, 3, 9368,
                                                                       2408, 9386, 98, 104, 2624,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9728, 0, 3, 9386,
                                                                       2417, 9404, 104, 110,
                                                                       2642, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9764, 0, 3, 9404,
                                                                       2426, 9422, 110, 116,
                                                                       2660, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9800, 0, 3, 9422,
                                                                       2435, 9440, 116, 122,
                                                                       2678, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9836, 0, 3, 9440,
                                                                       2444, 9458, 122, 128,
                                                                       2696, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9872, 0, 3, 9458,
                                                                       2453, 9476, 128, 134,
                                                                       2714, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9908, 0, 3, 9476,
                                                                       2462, 9494, 134, 140,
                                                                       2732, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9944, 0, 3, 9494,
                                                                       2471, 9512, 140, 146,
                                                                       2750, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9980, 0, 3, 9530,
                                                                       2507, 9548, 158, 164,
                                                                       2804, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10016, 0, 3, 9548,
                                                                       2516, 9566, 164, 170,
                                                                       2822, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10052, 0, 3, 9566,
                                                                       2525, 9584, 170, 176,
                                                                       2840, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10088, 0, 3, 9584,
                                                                       2534, 9602, 176, 182,
                                                                       2858, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10124, 0, 3, 9602,
                                                                       2543, 9620, 182, 188,
                                                                       2876, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10160, 0, 3, 9620,
                                                                       2552, 9638, 188, 194,
                                                                       2894, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10196, 0, 3, 9638,
                                                                       2561, 9656, 194, 200,
                                                                       2912, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10232, 0, 3, 9656,
                                                                       2570, 9674, 200, 206,
                                                                       2930, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10268, 0, 3, 9692,
                                                                       2624, 9728, 218, 228,
                                                                       3008, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10328, 0, 3, 9728,
                                                                       2642, 9764, 228, 238,
                                                                       3038, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10388, 0, 3, 9764,
                                                                       2660, 9800, 238, 248,
                                                                       3068, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10448, 0, 3, 9800,
                                                                       2678, 9836, 248, 258,
                                                                       3098, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10508, 0, 3, 9836,
                                                                       2696, 9872, 258, 268,
                                                                       3128, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10568, 0, 3, 9872,
                                                                       2714, 9908, 268, 278,
                                                                       3158, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10628, 0, 3, 9908,
                                                                       2732, 9944, 278, 288,
                                                                       3188, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10688, 0, 3, 9980,
                                                                       2804, 10016, 308, 318,
                                                                       3278, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10748, 0, 3,
                                                                       10016, 2822, 10052, 318,
                                                                       328, 3308, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10808, 0, 3,
                                                                       10052, 2840, 10088, 328,
                                                                       338, 3338, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10868, 0, 3,
                                                                       10088, 2858, 10124, 338,
                                                                       348, 3368, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10928, 0, 3,
                                                                       10124, 2876, 10160, 348,
                                                                       358, 3398, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10988, 0, 3,
                                                                       10160, 2894, 10196, 358,
                                                                       368, 3428, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11048, 0, 3,
                                                                       10196, 2912, 10232, 368,
                                                                       378, 3458, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11108, 0, 3,
                                                                       10268, 3008, 10328, 398,
                                                                       413, 3578, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11198, 0, 3,
                                                                       10328, 3038, 10388, 413,
                                                                       428, 3623, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11288, 0, 3,
                                                                       10388, 3068, 10448, 428,
                                                                       443, 3668, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11378, 0, 3,
                                                                       10448, 3098, 10508, 443,
                                                                       458, 3713, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11468, 0, 3,
                                                                       10508, 3128, 10568, 458,
                                                                       473, 3758, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11558, 0, 3,
                                                                       10568, 3158, 10628, 473,
                                                                       488, 3803, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11648, 0, 3,
                                                                       10688, 3278, 10748, 518,
                                                                       533, 3938, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11738, 0, 3,
                                                                       10748, 3308, 10808, 533,
                                                                       548, 3983, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11828, 0, 3,
                                                                       10808, 3338, 10868, 548,
                                                                       563, 4028, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11918, 0, 3,
                                                                       10868, 3368, 10928, 563,
                                                                       578, 4073, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12008, 0, 3,
                                                                       10928, 3398, 10988, 578,
                                                                       593, 4118, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12098, 0, 3,
                                                                       10988, 3428, 11048, 593,
                                                                       608, 4163, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12188, 0, 3,
                                                                       11108, 3578, 11198, 638,
                                                                       659, 4334, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12314, 0, 3,
                                                                       11198, 3623, 11288, 659,
                                                                       680, 4397, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12440, 0, 3,
                                                                       11288, 3668, 11378, 680,
                                                                       701, 4460, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12566, 0, 3,
                                                                       11378, 3713, 11468, 701,
                                                                       722, 4523, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12692, 0, 3,
                                                                       11468, 3758, 11558, 722,
                                                                       743, 4586, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12818, 0, 3,
                                                                       11648, 3938, 11738, 785,
                                                                       806, 4775, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12944, 0, 3,
                                                                       11738, 3983, 11828, 806,
                                                                       827, 4838, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13070, 0, 3,
                                                                       11828, 4028, 11918, 827,
                                                                       848, 4901, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13196, 0, 3,
                                                                       11918, 4073, 12008, 848,
                                                                       869, 4964, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13322, 0, 3,
                                                                       12008, 4118, 12098, 869,
                                                                       890, 5027, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13448, 0, 3,
                                                                       12188, 4334, 12314, 932,
                                                                       960, 5258, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13616, 0, 3,
                                                                       12314, 4397, 12440, 960,
                                                                       988, 5342, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13784, 0, 3,
                                                                       12440, 4460, 12566, 988,
                                                                       1016, 5426, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13952, 0, 3,
                                                                       12566, 4523, 12692, 1016,
                                                                       1044, 5510, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14120, 0, 3,
                                                                       12818, 4775, 12944, 1100,
                                                                       1128, 5762, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14288, 0, 3,
                                                                       12944, 4838, 13070, 1128,
                                                                       1156, 5846, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14456, 0, 3,
                                                                       13070, 4901, 13196, 1156,
                                                                       1184, 5930, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14624, 0, 3,
                                                                       13196, 4964, 13322, 1184,
                                                                       1212, 6014, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 14792, 0, 3,
                                                                       13448, 5258, 13616, 1268,
                                                                       1304, 6314, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15008, 0, 3,
                                                                       13616, 5342, 13784, 1304,
                                                                       1340, 6422, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15224, 0, 3,
                                                                       13784, 5426, 13952, 1340,
                                                                       1376, 6530, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15440, 0, 3,
                                                                       14120, 5762, 14288, 1448,
                                                                       1484, 6854, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15656, 0, 3,
                                                                       14288, 5846, 14456, 1484,
                                                                       1520, 6962, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15872, 0, 3,
                                                                       14456, 5930, 14624, 1520,
                                                                       1556, 7070, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 16088, 0, 3,
                                                                       14792, 6314, 15008, 1628,
                                                                       1673, 7448, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 16358, 0, 3,
                                                                       15008, 6422, 15224, 1673,
                                                                       1718, 7583, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 16628, 0, 3,
                                                                       15440, 6854, 15656, 1808,
                                                                       1853, 7988, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 16898, 0, 3,
                                                                       15656, 6962, 15872, 1853,
                                                                       1898, 8123, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 17168, 0, 3,
                                                                       16088, 7448, 16358, 1988,
                                                                       2043, 8588, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 17498, 0, 3,
                                                                       16628, 7988, 16898, 2153,
                                                                       2208, 9083, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17828, 3, 2318,
                                                                       2321, 9248, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17838, 3, 2321,
                                                                       2324, 9254, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17848, 3, 2324,
                                                                       2327, 9260, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17858, 3, 2327,
                                                                       2330, 9266, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17868, 3, 2330,
                                                                       2333, 9272, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17878, 3, 2333,
                                                                       2336, 9278, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17888, 3, 2336,
                                                                       2339, 9284, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17898, 3, 2339,
                                                                       2342, 9290, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17908, 3, 2342,
                                                                       2345, 9296, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17918, 3, 2345,
                                                                       2348, 9302, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17928, 3, 2354,
                                                                       2357, 9308, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17938, 3, 2357,
                                                                       2360, 9314, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17948, 3, 2360,
                                                                       2363, 9320, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17958, 3, 2363,
                                                                       2366, 9326, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17968, 3, 2366,
                                                                       2369, 9332, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17978, 3, 2369,
                                                                       2372, 9338, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17988, 3, 2372,
                                                                       2375, 9344, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17998, 3, 2375,
                                                                       2378, 9350, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18008, 3, 2378,
                                                                       2381, 9356, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18018, 3, 2381,
                                                                       2384, 9362, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18028, 0, 3,
                                                                       17828, 9248, 17838, 2390,
                                                                       2399, 9368, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18058, 0, 3,
                                                                       17838, 9254, 17848, 2399,
                                                                       2408, 9386, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18088, 0, 3,
                                                                       17848, 9260, 17858, 2408,
                                                                       2417, 9404, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18118, 0, 3,
                                                                       17858, 9266, 17868, 2417,
                                                                       2426, 9422, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18148, 0, 3,
                                                                       17868, 9272, 17878, 2426,
                                                                       2435, 9440, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18178, 0, 3,
                                                                       17878, 9278, 17888, 2435,
                                                                       2444, 9458, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18208, 0, 3,
                                                                       17888, 9284, 17898, 2444,
                                                                       2453, 9476, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18238, 0, 3,
                                                                       17898, 9290, 17908, 2453,
                                                                       2462, 9494, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18268, 0, 3,
                                                                       17908, 9296, 17918, 2462,
                                                                       2471, 9512, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18298, 0, 3,
                                                                       17928, 9308, 17938, 2489,
                                                                       2498, 9530, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18328, 0, 3,
                                                                       17938, 9314, 17948, 2498,
                                                                       2507, 9548, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18358, 0, 3,
                                                                       17948, 9320, 17958, 2507,
                                                                       2516, 9566, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18388, 0, 3,
                                                                       17958, 9326, 17968, 2516,
                                                                       2525, 9584, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18418, 0, 3,
                                                                       17968, 9332, 17978, 2525,
                                                                       2534, 9602, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18448, 0, 3,
                                                                       17978, 9338, 17988, 2534,
                                                                       2543, 9620, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18478, 0, 3,
                                                                       17988, 9344, 17998, 2543,
                                                                       2552, 9638, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18508, 0, 3,
                                                                       17998, 9350, 18008, 2552,
                                                                       2561, 9656, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18538, 0, 3,
                                                                       18008, 9356, 18018, 2561,
                                                                       2570, 9674, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18568, 0, 3,
                                                                       18028, 9368, 18058, 2588,
                                                                       2606, 9692, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18628, 0, 3,
                                                                       18058, 9386, 18088, 2606,
                                                                       2624, 9728, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18688, 0, 3,
                                                                       18088, 9404, 18118, 2624,
                                                                       2642, 9764, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18748, 0, 3,
                                                                       18118, 9422, 18148, 2642,
                                                                       2660, 9800, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18808, 0, 3,
                                                                       18148, 9440, 18178, 2660,
                                                                       2678, 9836, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18868, 0, 3,
                                                                       18178, 9458, 18208, 2678,
                                                                       2696, 9872, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18928, 0, 3,
                                                                       18208, 9476, 18238, 2696,
                                                                       2714, 9908, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18988, 0, 3,
                                                                       18238, 9494, 18268, 2714,
                                                                       2732, 9944, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19048, 0, 3,
                                                                       18298, 9530, 18328, 2768,
                                                                       2786, 9980, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19108, 0, 3,
                                                                       18328, 9548, 18358, 2786,
                                                                       2804, 10016, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19168, 0, 3,
                                                                       18358, 9566, 18388, 2804,
                                                                       2822, 10052, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19228, 0, 3,
                                                                       18388, 9584, 18418, 2822,
                                                                       2840, 10088, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19288, 0, 3,
                                                                       18418, 9602, 18448, 2840,
                                                                       2858, 10124, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19348, 0, 3,
                                                                       18448, 9620, 18478, 2858,
                                                                       2876, 10160, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19408, 0, 3,
                                                                       18478, 9638, 18508, 2876,
                                                                       2894, 10196, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19468, 0, 3,
                                                                       18508, 9656, 18538, 2894,
                                                                       2912, 10232, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19528, 0, 3,
                                                                       18568, 9692, 18628, 2948,
                                                                       2978, 10268, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19628, 0, 3,
                                                                       18628, 9728, 18688, 2978,
                                                                       3008, 10328, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19728, 0, 3,
                                                                       18688, 9764, 18748, 3008,
                                                                       3038, 10388, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19828, 0, 3,
                                                                       18748, 9800, 18808, 3038,
                                                                       3068, 10448, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19928, 0, 3,
                                                                       18808, 9836, 18868, 3068,
                                                                       3098, 10508, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 20028, 0, 3,
                                                                       18868, 9872, 18928, 3098,
                                                                       3128, 10568, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 20128, 0, 3,
                                                                       18928, 9908, 18988, 3128,
                                                                       3158, 10628, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 20228, 0, 3,
                                                                       19048, 9980, 19108, 3218,
                                                                       3248, 10688, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 20328, 0, 3,
                                                                       19108, 10016, 19168, 3248,
                                                                       3278, 10748, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 20428, 0, 3,
                                                                       19168, 10052, 19228, 3278,
                                                                       3308, 10808, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 20528, 0, 3,
                                                                       19228, 10088, 19288, 3308,
                                                                       3338, 10868, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 20628, 0, 3,
                                                                       19288, 10124, 19348, 3338,
                                                                       3368, 10928, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 20728, 0, 3,
                                                                       19348, 10160, 19408, 3368,
                                                                       3398, 10988, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 20828, 0, 3,
                                                                       19408, 10196, 19468, 3398,
                                                                       3428, 11048, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20928, 0, 3,
                                                                       19528, 10268, 19628, 3488,
                                                                       3533, 11108, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 21078, 0, 3,
                                                                       19628, 10328, 19728, 3533,
                                                                       3578, 11198, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 21228, 0, 3,
                                                                       19728, 10388, 19828, 3578,
                                                                       3623, 11288, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 21378, 0, 3,
                                                                       19828, 10448, 19928, 3623,
                                                                       3668, 11378, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 21528, 0, 3,
                                                                       19928, 10508, 20028, 3668,
                                                                       3713, 11468, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 21678, 0, 3,
                                                                       20028, 10568, 20128, 3713,
                                                                       3758, 11558, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 21828, 0, 3,
                                                                       20228, 10688, 20328, 3848,
                                                                       3893, 11648, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 21978, 0, 3,
                                                                       20328, 10748, 20428, 3893,
                                                                       3938, 11738, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 22128, 0, 3,
                                                                       20428, 10808, 20528, 3938,
                                                                       3983, 11828, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 22278, 0, 3,
                                                                       20528, 10868, 20628, 3983,
                                                                       4028, 11918, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 22428, 0, 3,
                                                                       20628, 10928, 20728, 4028,
                                                                       4073, 12008, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 22578, 0, 3,
                                                                       20728, 10988, 20828, 4073,
                                                                       4118, 12098, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 22728, 0, 3,
                                                                       20928, 11108, 21078, 4208,
                                                                       4271, 12188, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 22938, 0, 3,
                                                                       21078, 11198, 21228, 4271,
                                                                       4334, 12314, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 23148, 0, 3,
                                                                       21228, 11288, 21378, 4334,
                                                                       4397, 12440, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 23358, 0, 3,
                                                                       21378, 11378, 21528, 4397,
                                                                       4460, 12566, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 23568, 0, 3,
                                                                       21528, 11468, 21678, 4460,
                                                                       4523, 12692, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 23778, 0, 3,
                                                                       21828, 11648, 21978, 4649,
                                                                       4712, 12818, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 23988, 0, 3,
                                                                       21978, 11738, 22128, 4712,
                                                                       4775, 12944, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 24198, 0, 3,
                                                                       22128, 11828, 22278, 4775,
                                                                       4838, 13070, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 24408, 0, 3,
                                                                       22278, 11918, 22428, 4838,
                                                                       4901, 13196, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 24618, 0, 3,
                                                                       22428, 12008, 22578, 4901,
                                                                       4964, 13322, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 24828, 0, 3,
                                                                       22728, 12188, 22938, 5090,
                                                                       5174, 13448, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 25108, 0, 3,
                                                                       22938, 12314, 23148, 5174,
                                                                       5258, 13616, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 25388, 0, 3,
                                                                       23148, 12440, 23358, 5258,
                                                                       5342, 13784, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 25668, 0, 3,
                                                                       23358, 12566, 23568, 5342,
                                                                       5426, 13952, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 25948, 0, 3,
                                                                       23778, 12818, 23988, 5594,
                                                                       5678, 14120, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 26228, 0, 3,
                                                                       23988, 12944, 24198, 5678,
                                                                       5762, 14288, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 26508, 0, 3,
                                                                       24198, 13070, 24408, 5762,
                                                                       5846, 14456, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 26788, 0, 3,
                                                                       24408, 13196, 24618, 5846,
                                                                       5930, 14624, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 27068, 0, 3,
                                                                       24828, 13448, 25108, 6098,
                                                                       6206, 14792, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 27428, 0, 3,
                                                                       25108, 13616, 25388, 6206,
                                                                       6314, 15008, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 27788, 0, 3,
                                                                       25388, 13784, 25668, 6314,
                                                                       6422, 15224, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 28148, 0, 3,
                                                                       25948, 14120, 26228, 6638,
                                                                       6746, 15440, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 28508, 0, 3,
                                                                       26228, 14288, 26508, 6746,
                                                                       6854, 15656, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 28868, 0, 3,
                                                                       26508, 14456, 26788, 6854,
                                                                       6962, 15872, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 29228, 0, 3,
                                                                       27068, 14792, 27428, 7178,
                                                                       7313, 16088, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 29678, 0, 3,
                                                                       27428, 15008, 27788, 7313,
                                                                       7448, 16358, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 30128, 0, 3,
                                                                       28148, 15440, 28508, 7718,
                                                                       7853, 16628, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 30578, 0, 3,
                                                                       28508, 15656, 28868, 7853,
                                                                       7988, 16898, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 31028, 0, 3,
                                                                       29228, 16088, 29678, 8258,
                                                                       8423, 17168, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 31578, 0, 3,
                                                                       30128, 16628, 30578, 8753,
                                                                       8918, 17498, ncols, gamma,
                                                                       p, q);

                    simdfunc::contract_primitives(buffer, 32128, 24828, 280, ncols);

                    simdfunc::contract_primitives(buffer, 32604, 25948, 280, ncols);

                    simdfunc::contract_primitives(buffer, 33080, 27068, 360, ncols);

                    simdfunc::contract_primitives(buffer, 33692, 28148, 360, ncols);

                    simdfunc::contract_primitives(buffer, 34304, 29228, 450, ncols);

                    simdfunc::contract_primitives(buffer, 35069, 30128, 450, ncols);

                    simdfunc::contract_primitives(buffer, 35834, 31028, 550, ncols);

                    simdfunc::contract_primitives(buffer, 36769, 31578, 550, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 32408, 32128, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 32884, 32604, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 33440, 33080, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 34052, 33692, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 34754, 34304, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 35519, 35069, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 36384, 35834, 55, 1, nmax);

        simdtrf::transform_f_inner(buffer, 37319, 36769, 55, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 37704, 32408, 33440, 7, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 38292, 32884, 34052, 7, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 38880, 33440, 34754, 7, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 39636, 34052, 35519, 7, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 40392, 34754, 36384, 7, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 41337, 35519, 37319, 7, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 42282, 37704, 38880, 7, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 43458, 38292, 39636, 7, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 44634, 38880, 40392, 7, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 46146, 39636, 41337, 7, nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 47658, 42282, 44634, 7, nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 49618, 43458, 46146, 7, nmax);

        simdtrf::transform_f_inner(buffer, 51578, 49618, 28, 7, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 51578, 49, nmax);

        simdtrf::transform_f_inner(buffer, 51578, 47658, 28, 7, nmax);

        simdtrf::transform_i_outer(values + 637 * nvalues + n * npairs, nvalues, buffer, 51578,
                                   49, nmax);
    }

    for (size_t m = 0; m < 1274; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
