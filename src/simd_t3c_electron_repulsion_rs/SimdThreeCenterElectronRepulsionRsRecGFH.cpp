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


#include "SimdThreeCenterElectronRepulsionRsRecGFH.hpp"

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
#include "SimdTransferGD.hpp"
#include "SimdTransferGF.hpp"
#include "SimdTransferGP.hpp"
#include "SimdTransferHD.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_gfh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_gfh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 68319, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1386 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 68319, 48488, 6004, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1628, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1631, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1634, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1637, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1640, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1643, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1646, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1649, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1652, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1655, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1658, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1661, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1664, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1667, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1670, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1673, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1676, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1679, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1682, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1685, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1688, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1691, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1694, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1697, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1700, 3, 7, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1709, 3, 8, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1718, 3, 9, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1727, 3, 10, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1736, 3, 11, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1745, 3, 12, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1754, 3, 13, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1763, 3, 14, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1772, 3, 15, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1781, 3, 16, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1790, 3, 17, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1799, 3, 20, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1808, 3, 21, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1817, 3, 22, 71,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1826, 3, 23, 74,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1835, 3, 24, 77,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1844, 3, 25, 80,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1853, 3, 26, 83,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1862, 3, 27, 86,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1871, 3, 28, 89,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1880, 3, 29, 92,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1889, 3, 30, 95,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1898, 3, 32, 98,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1916, 3, 35, 104,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1934, 3, 38, 110,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1952, 3, 41, 116,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1970, 3, 44, 122,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1988, 3, 47, 128,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2006, 3, 50, 134,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2024, 3, 53, 140,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2042, 3, 56, 146,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2060, 3, 59, 152,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2078, 3, 65, 158,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2096, 3, 68, 164,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2114, 3, 71, 170,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2132, 3, 74, 176,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2150, 3, 77, 182,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2168, 3, 80, 188,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2186, 3, 83, 194,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2204, 3, 86, 200,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2222, 3, 89, 206,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2240, 3, 92, 212,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2258, 3, 98, 218,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2288, 3, 104, 228,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2318, 3, 110, 238,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2348, 3, 116, 248,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2378, 3, 122, 258,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2408, 3, 128, 268,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2438, 3, 134, 278,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2468, 3, 140, 288,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2498, 3, 146, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2528, 3, 158, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2558, 3, 164, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2588, 3, 170, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2618, 3, 176, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2648, 3, 182, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2678, 3, 188, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2708, 3, 194, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2738, 3, 200, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2768, 3, 206, 388,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2798, 3, 218, 398,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2843, 3, 228, 413,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2888, 3, 238, 428,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2933, 3, 248, 443,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2978, 3, 258, 458,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3023, 3, 268, 473,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3068, 3, 278, 488,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3113, 3, 288, 503,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3158, 3, 308, 518,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3203, 3, 318, 533,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3248, 3, 328, 548,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3293, 3, 338, 563,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3338, 3, 348, 578,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3383, 3, 358, 593,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3428, 3, 368, 608,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3473, 3, 378, 623,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3518, 3, 398, 638,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3581, 3, 413, 659,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3644, 3, 428, 680,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3707, 3, 443, 701,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3770, 3, 458, 722,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3833, 3, 473, 743,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3896, 3, 488, 764,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3959, 3, 518, 785,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4022, 3, 533, 806,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4085, 3, 548, 827,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4148, 3, 563, 848,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4211, 3, 578, 869,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4274, 3, 593, 890,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4337, 3, 608, 911,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4400, 3, 638, 932,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4484, 3, 659, 960,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4568, 3, 680, 988,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4652, 3, 701,
                                                                       1016, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4736, 3, 722,
                                                                       1044, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4820, 3, 743,
                                                                       1072, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4904, 3, 785,
                                                                       1100, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4988, 3, 806,
                                                                       1128, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5072, 3, 827,
                                                                       1156, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5156, 3, 848,
                                                                       1184, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5240, 3, 869,
                                                                       1212, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5324, 3, 890,
                                                                       1240, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5408, 3, 932,
                                                                       1268, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5516, 3, 960,
                                                                       1304, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5624, 3, 988,
                                                                       1340, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5732, 3, 1016,
                                                                       1376, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5840, 3, 1044,
                                                                       1412, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5948, 3, 1100,
                                                                       1448, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6056, 3, 1128,
                                                                       1484, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6164, 3, 1156,
                                                                       1520, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6272, 3, 1184,
                                                                       1556, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6380, 3, 1212,
                                                                       1592, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6488, 3, 7, 8,
                                                                       1634, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6494, 3, 8, 9,
                                                                       1637, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6500, 3, 9, 10,
                                                                       1640, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6506, 3, 10, 11,
                                                                       1643, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6512, 3, 11, 12,
                                                                       1646, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6518, 3, 12, 13,
                                                                       1649, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6524, 3, 13, 14,
                                                                       1652, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6530, 3, 14, 15,
                                                                       1655, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6536, 3, 15, 16,
                                                                       1658, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6542, 3, 16, 17,
                                                                       1661, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6548, 3, 20, 21,
                                                                       1670, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6554, 3, 21, 22,
                                                                       1673, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6560, 3, 22, 23,
                                                                       1676, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6566, 3, 23, 24,
                                                                       1679, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6572, 3, 24, 25,
                                                                       1682, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6578, 3, 25, 26,
                                                                       1685, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6584, 3, 26, 27,
                                                                       1688, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6590, 3, 27, 28,
                                                                       1691, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6596, 3, 28, 29,
                                                                       1694, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6602, 3, 29, 30,
                                                                       1697, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6608, 0, 3, 6488,
                                                                       1634, 6494, 1718, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6626, 0, 3, 6494,
                                                                       1637, 6500, 1727, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6644, 0, 3, 6500,
                                                                       1640, 6506, 1736, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6662, 0, 3, 6506,
                                                                       1643, 6512, 1745, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6680, 0, 3, 6512,
                                                                       1646, 6518, 1754, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6698, 0, 3, 6518,
                                                                       1649, 6524, 1763, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6716, 0, 3, 6524,
                                                                       1652, 6530, 1772, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6734, 0, 3, 6530,
                                                                       1655, 6536, 1781, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6752, 0, 3, 6536,
                                                                       1658, 6542, 1790, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6770, 0, 3, 6548,
                                                                       1670, 6554, 1817, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6788, 0, 3, 6554,
                                                                       1673, 6560, 1826, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6806, 0, 3, 6560,
                                                                       1676, 6566, 1835, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6824, 0, 3, 6566,
                                                                       1679, 6572, 1844, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6842, 0, 3, 6572,
                                                                       1682, 6578, 1853, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6860, 0, 3, 6578,
                                                                       1685, 6584, 1862, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6878, 0, 3, 6584,
                                                                       1688, 6590, 1871, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6896, 0, 3, 6590,
                                                                       1691, 6596, 1880, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6914, 0, 3, 6596,
                                                                       1694, 6602, 1889, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6932, 0, 3, 6608,
                                                                       1718, 6626, 98, 104, 1934,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6968, 0, 3, 6626,
                                                                       1727, 6644, 104, 110,
                                                                       1952, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7004, 0, 3, 6644,
                                                                       1736, 6662, 110, 116,
                                                                       1970, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7040, 0, 3, 6662,
                                                                       1745, 6680, 116, 122,
                                                                       1988, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7076, 0, 3, 6680,
                                                                       1754, 6698, 122, 128,
                                                                       2006, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7112, 0, 3, 6698,
                                                                       1763, 6716, 128, 134,
                                                                       2024, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7148, 0, 3, 6716,
                                                                       1772, 6734, 134, 140,
                                                                       2042, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7184, 0, 3, 6734,
                                                                       1781, 6752, 140, 146,
                                                                       2060, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7220, 0, 3, 6770,
                                                                       1817, 6788, 158, 164,
                                                                       2114, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7256, 0, 3, 6788,
                                                                       1826, 6806, 164, 170,
                                                                       2132, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7292, 0, 3, 6806,
                                                                       1835, 6824, 170, 176,
                                                                       2150, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7328, 0, 3, 6824,
                                                                       1844, 6842, 176, 182,
                                                                       2168, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7364, 0, 3, 6842,
                                                                       1853, 6860, 182, 188,
                                                                       2186, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7400, 0, 3, 6860,
                                                                       1862, 6878, 188, 194,
                                                                       2204, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7436, 0, 3, 6878,
                                                                       1871, 6896, 194, 200,
                                                                       2222, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7472, 0, 3, 6896,
                                                                       1880, 6914, 200, 206,
                                                                       2240, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7508, 0, 3, 6932,
                                                                       1934, 6968, 218, 228,
                                                                       2318, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7568, 0, 3, 6968,
                                                                       1952, 7004, 228, 238,
                                                                       2348, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7628, 0, 3, 7004,
                                                                       1970, 7040, 238, 248,
                                                                       2378, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7688, 0, 3, 7040,
                                                                       1988, 7076, 248, 258,
                                                                       2408, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7748, 0, 3, 7076,
                                                                       2006, 7112, 258, 268,
                                                                       2438, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7808, 0, 3, 7112,
                                                                       2024, 7148, 268, 278,
                                                                       2468, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7868, 0, 3, 7148,
                                                                       2042, 7184, 278, 288,
                                                                       2498, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7928, 0, 3, 7220,
                                                                       2114, 7256, 308, 318,
                                                                       2588, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7988, 0, 3, 7256,
                                                                       2132, 7292, 318, 328,
                                                                       2618, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8048, 0, 3, 7292,
                                                                       2150, 7328, 328, 338,
                                                                       2648, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8108, 0, 3, 7328,
                                                                       2168, 7364, 338, 348,
                                                                       2678, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8168, 0, 3, 7364,
                                                                       2186, 7400, 348, 358,
                                                                       2708, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8228, 0, 3, 7400,
                                                                       2204, 7436, 358, 368,
                                                                       2738, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8288, 0, 3, 7436,
                                                                       2222, 7472, 368, 378,
                                                                       2768, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8348, 0, 3, 7508,
                                                                       2318, 7568, 398, 413,
                                                                       2888, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8438, 0, 3, 7568,
                                                                       2348, 7628, 413, 428,
                                                                       2933, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8528, 0, 3, 7628,
                                                                       2378, 7688, 428, 443,
                                                                       2978, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8618, 0, 3, 7688,
                                                                       2408, 7748, 443, 458,
                                                                       3023, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8708, 0, 3, 7748,
                                                                       2438, 7808, 458, 473,
                                                                       3068, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8798, 0, 3, 7808,
                                                                       2468, 7868, 473, 488,
                                                                       3113, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8888, 0, 3, 7928,
                                                                       2588, 7988, 518, 533,
                                                                       3248, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8978, 0, 3, 7988,
                                                                       2618, 8048, 533, 548,
                                                                       3293, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9068, 0, 3, 8048,
                                                                       2648, 8108, 548, 563,
                                                                       3338, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9158, 0, 3, 8108,
                                                                       2678, 8168, 563, 578,
                                                                       3383, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9248, 0, 3, 8168,
                                                                       2708, 8228, 578, 593,
                                                                       3428, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9338, 0, 3, 8228,
                                                                       2738, 8288, 593, 608,
                                                                       3473, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9428, 0, 3, 8348,
                                                                       2888, 8438, 638, 659,
                                                                       3644, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9554, 0, 3, 8438,
                                                                       2933, 8528, 659, 680,
                                                                       3707, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9680, 0, 3, 8528,
                                                                       2978, 8618, 680, 701,
                                                                       3770, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9806, 0, 3, 8618,
                                                                       3023, 8708, 701, 722,
                                                                       3833, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9932, 0, 3, 8708,
                                                                       3068, 8798, 722, 743,
                                                                       3896, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10058, 0, 3, 8888,
                                                                       3248, 8978, 785, 806,
                                                                       4085, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10184, 0, 3, 8978,
                                                                       3293, 9068, 806, 827,
                                                                       4148, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10310, 0, 3, 9068,
                                                                       3338, 9158, 827, 848,
                                                                       4211, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10436, 0, 3, 9158,
                                                                       3383, 9248, 848, 869,
                                                                       4274, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10562, 0, 3, 9248,
                                                                       3428, 9338, 869, 890,
                                                                       4337, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10688, 0, 3, 9428,
                                                                       3644, 9554, 932, 960,
                                                                       4568, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10856, 0, 3, 9554,
                                                                       3707, 9680, 960, 988,
                                                                       4652, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11024, 0, 3, 9680,
                                                                       3770, 9806, 988, 1016,
                                                                       4736, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11192, 0, 3, 9806,
                                                                       3833, 9932, 1016, 1044,
                                                                       4820, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11360, 0, 3,
                                                                       10058, 4085, 10184, 1100,
                                                                       1128, 5072, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11528, 0, 3,
                                                                       10184, 4148, 10310, 1128,
                                                                       1156, 5156, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11696, 0, 3,
                                                                       10310, 4211, 10436, 1156,
                                                                       1184, 5240, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11864, 0, 3,
                                                                       10436, 4274, 10562, 1184,
                                                                       1212, 5324, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12032, 0, 3,
                                                                       10688, 4568, 10856, 1268,
                                                                       1304, 5624, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12248, 0, 3,
                                                                       10856, 4652, 11024, 1304,
                                                                       1340, 5732, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12464, 0, 3,
                                                                       11024, 4736, 11192, 1340,
                                                                       1376, 5840, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12680, 0, 3,
                                                                       11360, 5072, 11528, 1448,
                                                                       1484, 6164, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12896, 0, 3,
                                                                       11528, 5156, 11696, 1484,
                                                                       1520, 6272, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 13112, 0, 3,
                                                                       11696, 5240, 11864, 1520,
                                                                       1556, 6380, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13328, 3, 1628,
                                                                       1631, 6488, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13338, 3, 1631,
                                                                       1634, 6494, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13348, 3, 1634,
                                                                       1637, 6500, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13358, 3, 1637,
                                                                       1640, 6506, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13368, 3, 1640,
                                                                       1643, 6512, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13378, 3, 1643,
                                                                       1646, 6518, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13388, 3, 1646,
                                                                       1649, 6524, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13398, 3, 1649,
                                                                       1652, 6530, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13408, 3, 1652,
                                                                       1655, 6536, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13418, 3, 1655,
                                                                       1658, 6542, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13428, 3, 1664,
                                                                       1667, 6548, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13438, 3, 1667,
                                                                       1670, 6554, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13448, 3, 1670,
                                                                       1673, 6560, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13458, 3, 1673,
                                                                       1676, 6566, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13468, 3, 1676,
                                                                       1679, 6572, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13478, 3, 1679,
                                                                       1682, 6578, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13488, 3, 1682,
                                                                       1685, 6584, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13498, 3, 1685,
                                                                       1688, 6590, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13508, 3, 1688,
                                                                       1691, 6596, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13518, 3, 1691,
                                                                       1694, 6602, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13528, 0, 3,
                                                                       13328, 6488, 13338, 1700,
                                                                       1709, 6608, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13558, 0, 3,
                                                                       13338, 6494, 13348, 1709,
                                                                       1718, 6626, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13588, 0, 3,
                                                                       13348, 6500, 13358, 1718,
                                                                       1727, 6644, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13618, 0, 3,
                                                                       13358, 6506, 13368, 1727,
                                                                       1736, 6662, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13648, 0, 3,
                                                                       13368, 6512, 13378, 1736,
                                                                       1745, 6680, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13678, 0, 3,
                                                                       13378, 6518, 13388, 1745,
                                                                       1754, 6698, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13708, 0, 3,
                                                                       13388, 6524, 13398, 1754,
                                                                       1763, 6716, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13738, 0, 3,
                                                                       13398, 6530, 13408, 1763,
                                                                       1772, 6734, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13768, 0, 3,
                                                                       13408, 6536, 13418, 1772,
                                                                       1781, 6752, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13798, 0, 3,
                                                                       13428, 6548, 13438, 1799,
                                                                       1808, 6770, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13828, 0, 3,
                                                                       13438, 6554, 13448, 1808,
                                                                       1817, 6788, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13858, 0, 3,
                                                                       13448, 6560, 13458, 1817,
                                                                       1826, 6806, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13888, 0, 3,
                                                                       13458, 6566, 13468, 1826,
                                                                       1835, 6824, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13918, 0, 3,
                                                                       13468, 6572, 13478, 1835,
                                                                       1844, 6842, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13948, 0, 3,
                                                                       13478, 6578, 13488, 1844,
                                                                       1853, 6860, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13978, 0, 3,
                                                                       13488, 6584, 13498, 1853,
                                                                       1862, 6878, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14008, 0, 3,
                                                                       13498, 6590, 13508, 1862,
                                                                       1871, 6896, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14038, 0, 3,
                                                                       13508, 6596, 13518, 1871,
                                                                       1880, 6914, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14068, 0, 3,
                                                                       13528, 6608, 13558, 1898,
                                                                       1916, 6932, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14128, 0, 3,
                                                                       13558, 6626, 13588, 1916,
                                                                       1934, 6968, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14188, 0, 3,
                                                                       13588, 6644, 13618, 1934,
                                                                       1952, 7004, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14248, 0, 3,
                                                                       13618, 6662, 13648, 1952,
                                                                       1970, 7040, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14308, 0, 3,
                                                                       13648, 6680, 13678, 1970,
                                                                       1988, 7076, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14368, 0, 3,
                                                                       13678, 6698, 13708, 1988,
                                                                       2006, 7112, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14428, 0, 3,
                                                                       13708, 6716, 13738, 2006,
                                                                       2024, 7148, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14488, 0, 3,
                                                                       13738, 6734, 13768, 2024,
                                                                       2042, 7184, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14548, 0, 3,
                                                                       13798, 6770, 13828, 2078,
                                                                       2096, 7220, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14608, 0, 3,
                                                                       13828, 6788, 13858, 2096,
                                                                       2114, 7256, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14668, 0, 3,
                                                                       13858, 6806, 13888, 2114,
                                                                       2132, 7292, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14728, 0, 3,
                                                                       13888, 6824, 13918, 2132,
                                                                       2150, 7328, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14788, 0, 3,
                                                                       13918, 6842, 13948, 2150,
                                                                       2168, 7364, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14848, 0, 3,
                                                                       13948, 6860, 13978, 2168,
                                                                       2186, 7400, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14908, 0, 3,
                                                                       13978, 6878, 14008, 2186,
                                                                       2204, 7436, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14968, 0, 3,
                                                                       14008, 6896, 14038, 2204,
                                                                       2222, 7472, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15028, 0, 3,
                                                                       14068, 6932, 14128, 2258,
                                                                       2288, 7508, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15128, 0, 3,
                                                                       14128, 6968, 14188, 2288,
                                                                       2318, 7568, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15228, 0, 3,
                                                                       14188, 7004, 14248, 2318,
                                                                       2348, 7628, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15328, 0, 3,
                                                                       14248, 7040, 14308, 2348,
                                                                       2378, 7688, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15428, 0, 3,
                                                                       14308, 7076, 14368, 2378,
                                                                       2408, 7748, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15528, 0, 3,
                                                                       14368, 7112, 14428, 2408,
                                                                       2438, 7808, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15628, 0, 3,
                                                                       14428, 7148, 14488, 2438,
                                                                       2468, 7868, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15728, 0, 3,
                                                                       14548, 7220, 14608, 2528,
                                                                       2558, 7928, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15828, 0, 3,
                                                                       14608, 7256, 14668, 2558,
                                                                       2588, 7988, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15928, 0, 3,
                                                                       14668, 7292, 14728, 2588,
                                                                       2618, 8048, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 16028, 0, 3,
                                                                       14728, 7328, 14788, 2618,
                                                                       2648, 8108, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 16128, 0, 3,
                                                                       14788, 7364, 14848, 2648,
                                                                       2678, 8168, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 16228, 0, 3,
                                                                       14848, 7400, 14908, 2678,
                                                                       2708, 8228, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 16328, 0, 3,
                                                                       14908, 7436, 14968, 2708,
                                                                       2738, 8288, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 16428, 0, 3,
                                                                       15028, 7508, 15128, 2798,
                                                                       2843, 8348, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 16578, 0, 3,
                                                                       15128, 7568, 15228, 2843,
                                                                       2888, 8438, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 16728, 0, 3,
                                                                       15228, 7628, 15328, 2888,
                                                                       2933, 8528, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 16878, 0, 3,
                                                                       15328, 7688, 15428, 2933,
                                                                       2978, 8618, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17028, 0, 3,
                                                                       15428, 7748, 15528, 2978,
                                                                       3023, 8708, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17178, 0, 3,
                                                                       15528, 7808, 15628, 3023,
                                                                       3068, 8798, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17328, 0, 3,
                                                                       15728, 7928, 15828, 3158,
                                                                       3203, 8888, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17478, 0, 3,
                                                                       15828, 7988, 15928, 3203,
                                                                       3248, 8978, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17628, 0, 3,
                                                                       15928, 8048, 16028, 3248,
                                                                       3293, 9068, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17778, 0, 3,
                                                                       16028, 8108, 16128, 3293,
                                                                       3338, 9158, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17928, 0, 3,
                                                                       16128, 8168, 16228, 3338,
                                                                       3383, 9248, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18078, 0, 3,
                                                                       16228, 8228, 16328, 3383,
                                                                       3428, 9338, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 18228, 0, 3,
                                                                       16428, 8348, 16578, 3518,
                                                                       3581, 9428, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 18438, 0, 3,
                                                                       16578, 8438, 16728, 3581,
                                                                       3644, 9554, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 18648, 0, 3,
                                                                       16728, 8528, 16878, 3644,
                                                                       3707, 9680, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 18858, 0, 3,
                                                                       16878, 8618, 17028, 3707,
                                                                       3770, 9806, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19068, 0, 3,
                                                                       17028, 8708, 17178, 3770,
                                                                       3833, 9932, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19278, 0, 3,
                                                                       17328, 8888, 17478, 3959,
                                                                       4022, 10058, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19488, 0, 3,
                                                                       17478, 8978, 17628, 4022,
                                                                       4085, 10184, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19698, 0, 3,
                                                                       17628, 9068, 17778, 4085,
                                                                       4148, 10310, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19908, 0, 3,
                                                                       17778, 9158, 17928, 4148,
                                                                       4211, 10436, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 20118, 0, 3,
                                                                       17928, 9248, 18078, 4211,
                                                                       4274, 10562, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 20328, 0, 3,
                                                                       18228, 9428, 18438, 4400,
                                                                       4484, 10688, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 20608, 0, 3,
                                                                       18438, 9554, 18648, 4484,
                                                                       4568, 10856, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 20888, 0, 3,
                                                                       18648, 9680, 18858, 4568,
                                                                       4652, 11024, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 21168, 0, 3,
                                                                       18858, 9806, 19068, 4652,
                                                                       4736, 11192, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 21448, 0, 3,
                                                                       19278, 10058, 19488, 4904,
                                                                       4988, 11360, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 21728, 0, 3,
                                                                       19488, 10184, 19698, 4988,
                                                                       5072, 11528, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 22008, 0, 3,
                                                                       19698, 10310, 19908, 5072,
                                                                       5156, 11696, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 22288, 0, 3,
                                                                       19908, 10436, 20118, 5156,
                                                                       5240, 11864, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 22568, 0, 3,
                                                                       20328, 10688, 20608, 5408,
                                                                       5516, 12032, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 22928, 0, 3,
                                                                       20608, 10856, 20888, 5516,
                                                                       5624, 12248, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 23288, 0, 3,
                                                                       20888, 11024, 21168, 5624,
                                                                       5732, 12464, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 23648, 0, 3,
                                                                       21448, 11360, 21728, 5948,
                                                                       6056, 12680, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 24008, 0, 3,
                                                                       21728, 11528, 22008, 6056,
                                                                       6164, 12896, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 24368, 0, 3,
                                                                       22008, 11696, 22288, 6164,
                                                                       6272, 13112, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24728, 3, 6488,
                                                                       6494, 13348, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24743, 3, 6494,
                                                                       6500, 13358, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24758, 3, 6500,
                                                                       6506, 13368, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24773, 3, 6506,
                                                                       6512, 13378, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24788, 3, 6512,
                                                                       6518, 13388, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24803, 3, 6518,
                                                                       6524, 13398, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24818, 3, 6524,
                                                                       6530, 13408, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24833, 3, 6530,
                                                                       6536, 13418, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24848, 3, 6548,
                                                                       6554, 13448, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24863, 3, 6554,
                                                                       6560, 13458, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24878, 3, 6560,
                                                                       6566, 13468, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24893, 3, 6566,
                                                                       6572, 13478, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24908, 3, 6572,
                                                                       6578, 13488, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24923, 3, 6578,
                                                                       6584, 13498, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24938, 3, 6584,
                                                                       6590, 13508, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24953, 3, 6590,
                                                                       6596, 13518, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24968, 0, 3,
                                                                       24728, 13348, 24743, 6608,
                                                                       6626, 13588, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 25013, 0, 3,
                                                                       24743, 13358, 24758, 6626,
                                                                       6644, 13618, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 25058, 0, 3,
                                                                       24758, 13368, 24773, 6644,
                                                                       6662, 13648, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 25103, 0, 3,
                                                                       24773, 13378, 24788, 6662,
                                                                       6680, 13678, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 25148, 0, 3,
                                                                       24788, 13388, 24803, 6680,
                                                                       6698, 13708, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 25193, 0, 3,
                                                                       24803, 13398, 24818, 6698,
                                                                       6716, 13738, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 25238, 0, 3,
                                                                       24818, 13408, 24833, 6716,
                                                                       6734, 13768, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 25283, 0, 3,
                                                                       24848, 13448, 24863, 6770,
                                                                       6788, 13858, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 25328, 0, 3,
                                                                       24863, 13458, 24878, 6788,
                                                                       6806, 13888, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 25373, 0, 3,
                                                                       24878, 13468, 24893, 6806,
                                                                       6824, 13918, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 25418, 0, 3,
                                                                       24893, 13478, 24908, 6824,
                                                                       6842, 13948, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 25463, 0, 3,
                                                                       24908, 13488, 24923, 6842,
                                                                       6860, 13978, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 25508, 0, 3,
                                                                       24923, 13498, 24938, 6860,
                                                                       6878, 14008, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 25553, 0, 3,
                                                                       24938, 13508, 24953, 6878,
                                                                       6896, 14038, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25598, 0, 3,
                                                                       24968, 13588, 25013, 6932,
                                                                       6968, 14188, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25688, 0, 3,
                                                                       25013, 13618, 25058, 6968,
                                                                       7004, 14248, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25778, 0, 3,
                                                                       25058, 13648, 25103, 7004,
                                                                       7040, 14308, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25868, 0, 3,
                                                                       25103, 13678, 25148, 7040,
                                                                       7076, 14368, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25958, 0, 3,
                                                                       25148, 13708, 25193, 7076,
                                                                       7112, 14428, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 26048, 0, 3,
                                                                       25193, 13738, 25238, 7112,
                                                                       7148, 14488, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 26138, 0, 3,
                                                                       25283, 13858, 25328, 7220,
                                                                       7256, 14668, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 26228, 0, 3,
                                                                       25328, 13888, 25373, 7256,
                                                                       7292, 14728, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 26318, 0, 3,
                                                                       25373, 13918, 25418, 7292,
                                                                       7328, 14788, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 26408, 0, 3,
                                                                       25418, 13948, 25463, 7328,
                                                                       7364, 14848, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 26498, 0, 3,
                                                                       25463, 13978, 25508, 7364,
                                                                       7400, 14908, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 26588, 0, 3,
                                                                       25508, 14008, 25553, 7400,
                                                                       7436, 14968, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26678, 0, 3,
                                                                       25598, 14188, 25688, 7508,
                                                                       7568, 15228, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26828, 0, 3,
                                                                       25688, 14248, 25778, 7568,
                                                                       7628, 15328, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26978, 0, 3,
                                                                       25778, 14308, 25868, 7628,
                                                                       7688, 15428, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 27128, 0, 3,
                                                                       25868, 14368, 25958, 7688,
                                                                       7748, 15528, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 27278, 0, 3,
                                                                       25958, 14428, 26048, 7748,
                                                                       7808, 15628, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 27428, 0, 3,
                                                                       26138, 14668, 26228, 7928,
                                                                       7988, 15928, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 27578, 0, 3,
                                                                       26228, 14728, 26318, 7988,
                                                                       8048, 16028, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 27728, 0, 3,
                                                                       26318, 14788, 26408, 8048,
                                                                       8108, 16128, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 27878, 0, 3,
                                                                       26408, 14848, 26498, 8108,
                                                                       8168, 16228, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 28028, 0, 3,
                                                                       26498, 14908, 26588, 8168,
                                                                       8228, 16328, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 28178, 0, 3,
                                                                       26678, 15228, 26828, 8348,
                                                                       8438, 16728, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 28403, 0, 3,
                                                                       26828, 15328, 26978, 8438,
                                                                       8528, 16878, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 28628, 0, 3,
                                                                       26978, 15428, 27128, 8528,
                                                                       8618, 17028, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 28853, 0, 3,
                                                                       27128, 15528, 27278, 8618,
                                                                       8708, 17178, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 29078, 0, 3,
                                                                       27428, 15928, 27578, 8888,
                                                                       8978, 17628, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 29303, 0, 3,
                                                                       27578, 16028, 27728, 8978,
                                                                       9068, 17778, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 29528, 0, 3,
                                                                       27728, 16128, 27878, 9068,
                                                                       9158, 17928, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 29753, 0, 3,
                                                                       27878, 16228, 28028, 9158,
                                                                       9248, 18078, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 29978, 0, 3,
                                                                       28178, 16728, 28403, 9428,
                                                                       9554, 18648, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 30293, 0, 3,
                                                                       28403, 16878, 28628, 9554,
                                                                       9680, 18858, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 30608, 0, 3,
                                                                       28628, 17028, 28853, 9680,
                                                                       9806, 19068, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 30923, 0, 3,
                                                                       29078, 17628, 29303,
                                                                       10058, 10184, 19698,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 31238, 0, 3,
                                                                       29303, 17778, 29528,
                                                                       10184, 10310, 19908,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 31553, 0, 3,
                                                                       29528, 17928, 29753,
                                                                       10310, 10436, 20118,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 31868, 0, 3,
                                                                       29978, 18648, 30293,
                                                                       10688, 10856, 20888,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 32288, 0, 3,
                                                                       30293, 18858, 30608,
                                                                       10856, 11024, 21168,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 32708, 0, 3,
                                                                       30923, 19698, 31238,
                                                                       11360, 11528, 22008,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 33128, 0, 3,
                                                                       31238, 19908, 31553,
                                                                       11528, 11696, 22288,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 33548, 0, 3,
                                                                       31868, 20888, 32288,
                                                                       12032, 12248, 23288,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 34088, 0, 3,
                                                                       32708, 22008, 33128,
                                                                       12680, 12896, 24368,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34628, 3, 13328,
                                                                       13338, 24728, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34649, 3, 13338,
                                                                       13348, 24743, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34670, 3, 13348,
                                                                       13358, 24758, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34691, 3, 13358,
                                                                       13368, 24773, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34712, 3, 13368,
                                                                       13378, 24788, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34733, 3, 13378,
                                                                       13388, 24803, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34754, 3, 13388,
                                                                       13398, 24818, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34775, 3, 13398,
                                                                       13408, 24833, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34796, 3, 13428,
                                                                       13438, 24848, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34817, 3, 13438,
                                                                       13448, 24863, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34838, 3, 13448,
                                                                       13458, 24878, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34859, 3, 13458,
                                                                       13468, 24893, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34880, 3, 13468,
                                                                       13478, 24908, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34901, 3, 13478,
                                                                       13488, 24923, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34922, 3, 13488,
                                                                       13498, 24938, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34943, 3, 13498,
                                                                       13508, 24953, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 34964, 0, 3,
                                                                       34628, 24728, 34649,
                                                                       13528, 13558, 24968,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35027, 0, 3,
                                                                       34649, 24743, 34670,
                                                                       13558, 13588, 25013,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35090, 0, 3,
                                                                       34670, 24758, 34691,
                                                                       13588, 13618, 25058,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35153, 0, 3,
                                                                       34691, 24773, 34712,
                                                                       13618, 13648, 25103,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35216, 0, 3,
                                                                       34712, 24788, 34733,
                                                                       13648, 13678, 25148,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35279, 0, 3,
                                                                       34733, 24803, 34754,
                                                                       13678, 13708, 25193,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35342, 0, 3,
                                                                       34754, 24818, 34775,
                                                                       13708, 13738, 25238,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35405, 0, 3,
                                                                       34796, 24848, 34817,
                                                                       13798, 13828, 25283,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35468, 0, 3,
                                                                       34817, 24863, 34838,
                                                                       13828, 13858, 25328,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35531, 0, 3,
                                                                       34838, 24878, 34859,
                                                                       13858, 13888, 25373,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35594, 0, 3,
                                                                       34859, 24893, 34880,
                                                                       13888, 13918, 25418,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35657, 0, 3,
                                                                       34880, 24908, 34901,
                                                                       13918, 13948, 25463,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35720, 0, 3,
                                                                       34901, 24923, 34922,
                                                                       13948, 13978, 25508,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 35783, 0, 3,
                                                                       34922, 24938, 34943,
                                                                       13978, 14008, 25553,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 35846, 0, 3,
                                                                       34964, 24968, 35027,
                                                                       14068, 14128, 25598,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 35972, 0, 3,
                                                                       35027, 25013, 35090,
                                                                       14128, 14188, 25688,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 36098, 0, 3,
                                                                       35090, 25058, 35153,
                                                                       14188, 14248, 25778,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 36224, 0, 3,
                                                                       35153, 25103, 35216,
                                                                       14248, 14308, 25868,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 36350, 0, 3,
                                                                       35216, 25148, 35279,
                                                                       14308, 14368, 25958,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 36476, 0, 3,
                                                                       35279, 25193, 35342,
                                                                       14368, 14428, 26048,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 36602, 0, 3,
                                                                       35405, 25283, 35468,
                                                                       14548, 14608, 26138,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 36728, 0, 3,
                                                                       35468, 25328, 35531,
                                                                       14608, 14668, 26228,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 36854, 0, 3,
                                                                       35531, 25373, 35594,
                                                                       14668, 14728, 26318,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 36980, 0, 3,
                                                                       35594, 25418, 35657,
                                                                       14728, 14788, 26408,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 37106, 0, 3,
                                                                       35657, 25463, 35720,
                                                                       14788, 14848, 26498,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 37232, 0, 3,
                                                                       35720, 25508, 35783,
                                                                       14848, 14908, 26588,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 37358, 0, 3,
                                                                       35846, 25598, 35972,
                                                                       15028, 15128, 26678,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 37568, 0, 3,
                                                                       35972, 25688, 36098,
                                                                       15128, 15228, 26828,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 37778, 0, 3,
                                                                       36098, 25778, 36224,
                                                                       15228, 15328, 26978,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 37988, 0, 3,
                                                                       36224, 25868, 36350,
                                                                       15328, 15428, 27128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 38198, 0, 3,
                                                                       36350, 25958, 36476,
                                                                       15428, 15528, 27278,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 38408, 0, 3,
                                                                       36602, 26138, 36728,
                                                                       15728, 15828, 27428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 38618, 0, 3,
                                                                       36728, 26228, 36854,
                                                                       15828, 15928, 27578,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 38828, 0, 3,
                                                                       36854, 26318, 36980,
                                                                       15928, 16028, 27728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 39038, 0, 3,
                                                                       36980, 26408, 37106,
                                                                       16028, 16128, 27878,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 39248, 0, 3,
                                                                       37106, 26498, 37232,
                                                                       16128, 16228, 28028,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 39458, 0, 3,
                                                                       37358, 26678, 37568,
                                                                       16428, 16578, 28178,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 39773, 0, 3,
                                                                       37568, 26828, 37778,
                                                                       16578, 16728, 28403,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 40088, 0, 3,
                                                                       37778, 26978, 37988,
                                                                       16728, 16878, 28628,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 40403, 0, 3,
                                                                       37988, 27128, 38198,
                                                                       16878, 17028, 28853,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 40718, 0, 3,
                                                                       38408, 27428, 38618,
                                                                       17328, 17478, 29078,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 41033, 0, 3,
                                                                       38618, 27578, 38828,
                                                                       17478, 17628, 29303,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 41348, 0, 3,
                                                                       38828, 27728, 39038,
                                                                       17628, 17778, 29528,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 41663, 0, 3,
                                                                       39038, 27878, 39248,
                                                                       17778, 17928, 29753,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 41978, 0, 3,
                                                                       39458, 28178, 39773,
                                                                       18228, 18438, 29978,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 42419, 0, 3,
                                                                       39773, 28403, 40088,
                                                                       18438, 18648, 30293,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 42860, 0, 3,
                                                                       40088, 28628, 40403,
                                                                       18648, 18858, 30608,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 43301, 0, 3,
                                                                       40718, 29078, 41033,
                                                                       19278, 19488, 30923,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 43742, 0, 3,
                                                                       41033, 29303, 41348,
                                                                       19488, 19698, 31238,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 44183, 0, 3,
                                                                       41348, 29528, 41663,
                                                                       19698, 19908, 31553,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 44624, 0, 3,
                                                                       41978, 29978, 42419,
                                                                       20328, 20608, 31868,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 45212, 0, 3,
                                                                       42419, 30293, 42860,
                                                                       20608, 20888, 32288,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 45800, 0, 3,
                                                                       43301, 30923, 43742,
                                                                       21448, 21728, 32708,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 46388, 0, 3,
                                                                       43742, 31238, 44183,
                                                                       21728, 22008, 33128,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 46976, 0, 3,
                                                                       44624, 31868, 45212,
                                                                       22568, 22928, 33548,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 47732, 0, 3,
                                                                       45800, 32708, 46388,
                                                                       23648, 24008, 34088,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 48488, 39458, 315, ncols);

                    simdfunc::contract_primitives(buffer, 48968, 40718, 315, ncols);

                    simdfunc::contract_primitives(buffer, 49448, 41978, 441, ncols);

                    simdfunc::contract_primitives(buffer, 50120, 43301, 441, ncols);

                    simdfunc::contract_primitives(buffer, 50792, 44624, 588, ncols);

                    simdfunc::contract_primitives(buffer, 51688, 45800, 588, ncols);

                    simdfunc::contract_primitives(buffer, 52584, 46976, 756, ncols);

                    simdfunc::contract_primitives(buffer, 53736, 47732, 756, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 48803, 48488, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 49283, 48968, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 49889, 49448, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 50561, 50120, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 51380, 50792, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 52276, 51688, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 53340, 52584, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 54492, 53736, 36, 1, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 54888, 48803, 49889, 11,
                                             nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 55383, 49283, 50561, 11,
                                             nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 55878, 49889, 51380, 11,
                                             nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 56571, 50561, 52276, 11,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 57264, 51380, 53340, 11,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 58188, 52276, 54492, 11,
                                             nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 59112, 54888, 55878, 11,
                                             nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 60102, 55383, 56571, 11,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 61092, 55878, 57264, 11,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 62478, 56571, 58188, 11,
                                             nmax);

        simdtrf::compute_hrr_gf_out_of_first(buffer, coordinates, 63864, 59112, 61092, 11,
                                             nmax);

        simdtrf::compute_hrr_gf_out_of_first(buffer, coordinates, 65514, 60102, 62478, 11,
                                             nmax);

        simdtrf::transform_f_inner(buffer, 67164, 65514, 15, 11, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 67164, 77, nmax);

        simdtrf::transform_f_inner(buffer, 67164, 63864, 15, 11, nmax);

        simdtrf::transform_g_outer(values + 693 * nvalues + n * npairs, nvalues, buffer, 67164,
                                   77, nmax);
    }

    for (size_t m = 0; m < 1386; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
