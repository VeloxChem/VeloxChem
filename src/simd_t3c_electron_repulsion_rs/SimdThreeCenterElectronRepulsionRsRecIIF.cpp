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


#include "SimdThreeCenterElectronRepulsionRsRecIIF.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
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
#include "SimdTransformF.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_iif_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_iif_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 204100, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2366 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 204100, 80088, 12929, dimensions);

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
                                                            15}, ncols, fj, i * nprim_b + j, fq,
                                                            omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 22, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
                                                        ncols, fj, i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 38, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 41, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 44, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 47, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 50, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 53, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 56, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 59, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 62, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 65, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 68, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 71, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 74, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 77, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 80, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 83, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 86, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 89, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 92, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 95, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 98, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 101, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 104, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 107, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 110, 0, 3, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 113, 0, 3, 34, 35,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 116, 0, 3, 35, 36,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 119, 0, 3, 36, 37,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 122, 0, 3, 7, 8,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 128, 0, 3, 8, 9,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 134, 0, 3, 9, 10,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 140, 0, 3, 10, 11,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 146, 0, 3, 11, 12,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 152, 0, 3, 12, 13,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 158, 0, 3, 13, 14,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 164, 0, 3, 14, 15,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 170, 0, 3, 15, 16,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 176, 0, 3, 16, 17,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 182, 0, 3, 17, 18,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 188, 0, 3, 18, 19,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 194, 0, 3, 19, 20,
                                                                       74, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 200, 0, 3, 23, 24,
                                                                       80, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 206, 0, 3, 24, 25,
                                                                       83, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 212, 0, 3, 25, 26,
                                                                       86, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 218, 0, 3, 26, 27,
                                                                       89, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 224, 0, 3, 27, 28,
                                                                       92, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 230, 0, 3, 28, 29,
                                                                       95, 98, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 236, 0, 3, 29, 30,
                                                                       98, 101, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 242, 0, 3, 30, 31,
                                                                       101, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 248, 0, 3, 31, 32,
                                                                       104, 107, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 254, 0, 3, 32, 33,
                                                                       107, 110, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 260, 0, 3, 33, 34,
                                                                       110, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 266, 0, 3, 34, 35,
                                                                       113, 116, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 272, 0, 3, 35, 36,
                                                                       116, 119, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 278, 0, 3, 38, 41,
                                                                       122, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 288, 0, 3, 41, 44,
                                                                       128, 134, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 298, 0, 3, 44, 47,
                                                                       134, 140, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 308, 0, 3, 47, 50,
                                                                       140, 146, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 318, 0, 3, 50, 53,
                                                                       146, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 328, 0, 3, 53, 56,
                                                                       152, 158, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 338, 0, 3, 56, 59,
                                                                       158, 164, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 348, 0, 3, 59, 62,
                                                                       164, 170, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 62, 65,
                                                                       170, 176, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 368, 0, 3, 65, 68,
                                                                       176, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 378, 0, 3, 68, 71,
                                                                       182, 188, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 388, 0, 3, 71, 74,
                                                                       188, 194, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 398, 0, 3, 80, 83,
                                                                       200, 206, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 408, 0, 3, 83, 86,
                                                                       206, 212, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 418, 0, 3, 86, 89,
                                                                       212, 218, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 89, 92,
                                                                       218, 224, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 438, 0, 3, 92, 95,
                                                                       224, 230, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 448, 0, 3, 95, 98,
                                                                       230, 236, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 458, 0, 3, 98,
                                                                       101, 236, 242, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 468, 0, 3, 101,
                                                                       104, 242, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 478, 0, 3, 104,
                                                                       107, 248, 254, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 488, 0, 3, 107,
                                                                       110, 254, 260, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 498, 0, 3, 110,
                                                                       113, 260, 266, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 508, 0, 3, 113,
                                                                       116, 266, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 518, 0, 3, 122,
                                                                       128, 278, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 533, 0, 3, 128,
                                                                       134, 288, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 548, 0, 3, 134,
                                                                       140, 298, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 563, 0, 3, 140,
                                                                       146, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 578, 0, 3, 146,
                                                                       152, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 593, 0, 3, 152,
                                                                       158, 328, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 608, 0, 3, 158,
                                                                       164, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 623, 0, 3, 164,
                                                                       170, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 638, 0, 3, 170,
                                                                       176, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 653, 0, 3, 176,
                                                                       182, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 668, 0, 3, 182,
                                                                       188, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 683, 0, 3, 200,
                                                                       206, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 698, 0, 3, 206,
                                                                       212, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 713, 0, 3, 212,
                                                                       218, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 728, 0, 3, 218,
                                                                       224, 428, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 743, 0, 3, 224,
                                                                       230, 438, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 758, 0, 3, 230,
                                                                       236, 448, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 773, 0, 3, 236,
                                                                       242, 458, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 788, 0, 3, 242,
                                                                       248, 468, 478, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 803, 0, 3, 248,
                                                                       254, 478, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 818, 0, 3, 254,
                                                                       260, 488, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 833, 0, 3, 260,
                                                                       266, 498, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 848, 0, 3, 278,
                                                                       288, 518, 533, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 869, 0, 3, 288,
                                                                       298, 533, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 890, 0, 3, 298,
                                                                       308, 548, 563, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 911, 0, 3, 308,
                                                                       318, 563, 578, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 932, 0, 3, 318,
                                                                       328, 578, 593, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 953, 0, 3, 328,
                                                                       338, 593, 608, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 974, 0, 3, 338,
                                                                       348, 608, 623, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 995, 0, 3, 348,
                                                                       358, 623, 638, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1016, 0, 3, 358,
                                                                       368, 638, 653, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1037, 0, 3, 368,
                                                                       378, 653, 668, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1058, 0, 3, 398,
                                                                       408, 683, 698, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1079, 0, 3, 408,
                                                                       418, 698, 713, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1100, 0, 3, 418,
                                                                       428, 713, 728, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1121, 0, 3, 428,
                                                                       438, 728, 743, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1142, 0, 3, 438,
                                                                       448, 743, 758, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1163, 0, 3, 448,
                                                                       458, 758, 773, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1184, 0, 3, 458,
                                                                       468, 773, 788, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1205, 0, 3, 468,
                                                                       478, 788, 803, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1226, 0, 3, 478,
                                                                       488, 803, 818, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1247, 0, 3, 488,
                                                                       498, 818, 833, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 518,
                                                                       533, 848, 869, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1296, 0, 3, 533,
                                                                       548, 869, 890, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1324, 0, 3, 548,
                                                                       563, 890, 911, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1352, 0, 3, 563,
                                                                       578, 911, 932, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 578,
                                                                       593, 932, 953, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1408, 0, 3, 593,
                                                                       608, 953, 974, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1436, 0, 3, 608,
                                                                       623, 974, 995, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1464, 0, 3, 623,
                                                                       638, 995, 1016, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 638,
                                                                       653, 1016, 1037, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 683,
                                                                       698, 1058, 1079, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1548, 0, 3, 698,
                                                                       713, 1079, 1100, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1576, 0, 3, 713,
                                                                       728, 1100, 1121, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1604, 0, 3, 728,
                                                                       743, 1121, 1142, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1632, 0, 3, 743,
                                                                       758, 1142, 1163, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1660, 0, 3, 758,
                                                                       773, 1163, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 773,
                                                                       788, 1184, 1205, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1716, 0, 3, 788,
                                                                       803, 1205, 1226, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1744, 0, 3, 803,
                                                                       818, 1226, 1247, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1772, 0, 3, 848,
                                                                       869, 1268, 1296, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1808, 0, 3, 869,
                                                                       890, 1296, 1324, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1844, 0, 3, 890,
                                                                       911, 1324, 1352, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1880, 0, 3, 911,
                                                                       932, 1352, 1380, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1916, 0, 3, 932,
                                                                       953, 1380, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1952, 0, 3, 953,
                                                                       974, 1408, 1436, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1988, 0, 3, 974,
                                                                       995, 1436, 1464, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2024, 0, 3, 995,
                                                                       1016, 1464, 1492, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2060, 0, 3, 1058,
                                                                       1079, 1520, 1548, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2096, 0, 3, 1079,
                                                                       1100, 1548, 1576, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2132, 0, 3, 1100,
                                                                       1121, 1576, 1604, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2168, 0, 3, 1121,
                                                                       1142, 1604, 1632, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2204, 0, 3, 1142,
                                                                       1163, 1632, 1660, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2240, 0, 3, 1163,
                                                                       1184, 1660, 1688, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2276, 0, 3, 1184,
                                                                       1205, 1688, 1716, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2312, 0, 3, 1205,
                                                                       1226, 1716, 1744, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2348, 0, 3, 1268,
                                                                       1296, 1772, 1808, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2393, 0, 3, 1296,
                                                                       1324, 1808, 1844, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2438, 0, 3, 1324,
                                                                       1352, 1844, 1880, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2483, 0, 3, 1352,
                                                                       1380, 1880, 1916, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2528, 0, 3, 1380,
                                                                       1408, 1916, 1952, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2573, 0, 3, 1408,
                                                                       1436, 1952, 1988, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2618, 0, 3, 1436,
                                                                       1464, 1988, 2024, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2663, 0, 3, 1520,
                                                                       1548, 2060, 2096, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2708, 0, 3, 1548,
                                                                       1576, 2096, 2132, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2753, 0, 3, 1576,
                                                                       1604, 2132, 2168, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2798, 0, 3, 1604,
                                                                       1632, 2168, 2204, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2843, 0, 3, 1632,
                                                                       1660, 2204, 2240, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2888, 0, 3, 1660,
                                                                       1688, 2240, 2276, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2933, 0, 3, 1688,
                                                                       1716, 2276, 2312, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2978, 0, 3, 1772,
                                                                       1808, 2348, 2393, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3033, 0, 3, 1808,
                                                                       1844, 2393, 2438, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3088, 0, 3, 1844,
                                                                       1880, 2438, 2483, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3143, 0, 3, 1880,
                                                                       1916, 2483, 2528, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3198, 0, 3, 1916,
                                                                       1952, 2528, 2573, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3253, 0, 3, 1952,
                                                                       1988, 2573, 2618, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3308, 0, 3, 2060,
                                                                       2096, 2663, 2708, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3363, 0, 3, 2096,
                                                                       2132, 2708, 2753, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3418, 0, 3, 2132,
                                                                       2168, 2753, 2798, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3473, 0, 3, 2168,
                                                                       2204, 2798, 2843, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3528, 0, 3, 2204,
                                                                       2240, 2843, 2888, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3583, 0, 3, 2240,
                                                                       2276, 2888, 2933, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3638, 0, 3, 2348,
                                                                       2393, 2978, 3033, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3704, 0, 3, 2393,
                                                                       2438, 3033, 3088, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3770, 0, 3, 2438,
                                                                       2483, 3088, 3143, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3836, 0, 3, 2483,
                                                                       2528, 3143, 3198, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3902, 0, 3, 2528,
                                                                       2573, 3198, 3253, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3968, 0, 3, 2663,
                                                                       2708, 3308, 3363, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4034, 0, 3, 2708,
                                                                       2753, 3363, 3418, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4100, 0, 3, 2753,
                                                                       2798, 3418, 3473, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4166, 0, 3, 2798,
                                                                       2843, 3473, 3528, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4232, 0, 3, 2843,
                                                                       2888, 3528, 3583, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4298, 0, 3, 2978,
                                                                       3033, 3638, 3704, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4376, 0, 3, 3033,
                                                                       3088, 3704, 3770, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4454, 0, 3, 3088,
                                                                       3143, 3770, 3836, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4532, 0, 3, 3143,
                                                                       3198, 3836, 3902, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4610, 0, 3, 3308,
                                                                       3363, 3968, 4034, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4688, 0, 3, 3363,
                                                                       3418, 4034, 4100, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4766, 0, 3, 3418,
                                                                       3473, 4100, 4166, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4844, 0, 3, 3473,
                                                                       3528, 4166, 4232, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4922, 0, 3, 3638,
                                                                       3704, 4298, 4376, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 5013, 0, 3, 3704,
                                                                       3770, 4376, 4454, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 5104, 0, 3, 3770,
                                                                       3836, 4454, 4532, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 5195, 0, 3, 3968,
                                                                       4034, 4610, 4688, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 5286, 0, 3, 4034,
                                                                       4100, 4688, 4766, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 5377, 0, 3, 4100,
                                                                       4166, 4766, 4844, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5468, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5471, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5474, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5477, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5480, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5483, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5486, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5489, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5492, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5495, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5498, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5501, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5504, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5507, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5510, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5513, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5516, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5519, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5522, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5525, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5528, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5531, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5534, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5537, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5540, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5543, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5546, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5549, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5552, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5555, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5558, 3, 7, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5567, 3, 8, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5576, 3, 9, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5585, 3, 10, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5594, 3, 11, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5603, 3, 12, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5612, 3, 13, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5621, 3, 14, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5630, 3, 15, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5639, 3, 16, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5648, 3, 17, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5657, 3, 18, 71,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5666, 3, 19, 74,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5675, 3, 20, 77,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5684, 3, 23, 80,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5693, 3, 24, 83,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5702, 3, 25, 86,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5711, 3, 26, 89,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5720, 3, 27, 92,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5729, 3, 28, 95,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5738, 3, 29, 98,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5747, 3, 30, 101,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5756, 3, 31, 104,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5765, 3, 32, 107,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5774, 3, 33, 110,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5783, 3, 34, 113,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5792, 3, 35, 116,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5801, 3, 36, 119,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5810, 3, 38, 122,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5828, 3, 41, 128,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5846, 3, 44, 134,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5864, 3, 47, 140,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5882, 3, 50, 146,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5900, 3, 53, 152,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5918, 3, 56, 158,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5936, 3, 59, 164,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5954, 3, 62, 170,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5972, 3, 65, 176,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5990, 3, 68, 182,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6008, 3, 71, 188,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6026, 3, 74, 194,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6044, 3, 80, 200,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6062, 3, 83, 206,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6080, 3, 86, 212,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6098, 3, 89, 218,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6116, 3, 92, 224,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6134, 3, 95, 230,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6152, 3, 98, 236,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6170, 3, 101, 242,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6188, 3, 104, 248,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6206, 3, 107, 254,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6224, 3, 110, 260,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6242, 3, 113, 266,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6260, 3, 116, 272,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6278, 3, 122, 278,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6308, 3, 128, 288,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6338, 3, 134, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6368, 3, 140, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6398, 3, 146, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6428, 3, 152, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6458, 3, 158, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6488, 3, 164, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6518, 3, 170, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6548, 3, 176, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6578, 3, 182, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6608, 3, 188, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6638, 3, 200, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6668, 3, 206, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6698, 3, 212, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6728, 3, 218, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6758, 3, 224, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6788, 3, 230, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6818, 3, 236, 458,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6848, 3, 242, 468,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6878, 3, 248, 478,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6908, 3, 254, 488,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6938, 3, 260, 498,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6968, 3, 266, 508,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6998, 3, 278, 518,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7043, 3, 288, 533,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7088, 3, 298, 548,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7133, 3, 308, 563,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7178, 3, 318, 578,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7223, 3, 328, 593,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7268, 3, 338, 608,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7313, 3, 348, 623,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7358, 3, 358, 638,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7403, 3, 368, 653,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7448, 3, 378, 668,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7493, 3, 398, 683,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7538, 3, 408, 698,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7583, 3, 418, 713,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7628, 3, 428, 728,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7673, 3, 438, 743,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7718, 3, 448, 758,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7763, 3, 458, 773,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7808, 3, 468, 788,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7853, 3, 478, 803,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7898, 3, 488, 818,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7943, 3, 498, 833,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7988, 3, 518, 848,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8051, 3, 533, 869,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8114, 3, 548, 890,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8177, 3, 563, 911,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8240, 3, 578, 932,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8303, 3, 593, 953,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8366, 3, 608, 974,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8429, 3, 623, 995,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8492, 3, 638,
                                                                       1016, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8555, 3, 653,
                                                                       1037, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8618, 3, 683,
                                                                       1058, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8681, 3, 698,
                                                                       1079, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8744, 3, 713,
                                                                       1100, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8807, 3, 728,
                                                                       1121, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8870, 3, 743,
                                                                       1142, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8933, 3, 758,
                                                                       1163, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8996, 3, 773,
                                                                       1184, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9059, 3, 788,
                                                                       1205, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9122, 3, 803,
                                                                       1226, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 9185, 3, 818,
                                                                       1247, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9248, 3, 848,
                                                                       1268, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9332, 3, 869,
                                                                       1296, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9416, 3, 890,
                                                                       1324, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9500, 3, 911,
                                                                       1352, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9584, 3, 932,
                                                                       1380, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9668, 3, 953,
                                                                       1408, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9752, 3, 974,
                                                                       1436, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9836, 3, 995,
                                                                       1464, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9920, 3, 1016,
                                                                       1492, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10004, 3, 1058,
                                                                       1520, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10088, 3, 1079,
                                                                       1548, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10172, 3, 1100,
                                                                       1576, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10256, 3, 1121,
                                                                       1604, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10340, 3, 1142,
                                                                       1632, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10424, 3, 1163,
                                                                       1660, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10508, 3, 1184,
                                                                       1688, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10592, 3, 1205,
                                                                       1716, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 10676, 3, 1226,
                                                                       1744, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10760, 3, 1268,
                                                                       1772, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10868, 3, 1296,
                                                                       1808, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10976, 3, 1324,
                                                                       1844, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11084, 3, 1352,
                                                                       1880, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11192, 3, 1380,
                                                                       1916, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11300, 3, 1408,
                                                                       1952, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11408, 3, 1436,
                                                                       1988, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11516, 3, 1464,
                                                                       2024, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11624, 3, 1520,
                                                                       2060, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11732, 3, 1548,
                                                                       2096, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11840, 3, 1576,
                                                                       2132, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11948, 3, 1604,
                                                                       2168, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12056, 3, 1632,
                                                                       2204, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12164, 3, 1660,
                                                                       2240, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12272, 3, 1688,
                                                                       2276, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 12380, 3, 1716,
                                                                       2312, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12488, 3, 1772,
                                                                       2348, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12623, 3, 1808,
                                                                       2393, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12758, 3, 1844,
                                                                       2438, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12893, 3, 1880,
                                                                       2483, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13028, 3, 1916,
                                                                       2528, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13163, 3, 1952,
                                                                       2573, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13298, 3, 1988,
                                                                       2618, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13433, 3, 2060,
                                                                       2663, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13568, 3, 2096,
                                                                       2708, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13703, 3, 2132,
                                                                       2753, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13838, 3, 2168,
                                                                       2798, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13973, 3, 2204,
                                                                       2843, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14108, 3, 2240,
                                                                       2888, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 14243, 3, 2276,
                                                                       2933, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14378, 3, 2348,
                                                                       2978, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14543, 3, 2393,
                                                                       3033, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14708, 3, 2438,
                                                                       3088, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14873, 3, 2483,
                                                                       3143, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15038, 3, 2528,
                                                                       3198, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15203, 3, 2573,
                                                                       3253, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15368, 3, 2663,
                                                                       3308, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15533, 3, 2708,
                                                                       3363, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15698, 3, 2753,
                                                                       3418, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15863, 3, 2798,
                                                                       3473, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16028, 3, 2843,
                                                                       3528, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 16193, 3, 2888,
                                                                       3583, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 16358, 3, 2978,
                                                                       3638, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 16556, 3, 3033,
                                                                       3704, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 16754, 3, 3088,
                                                                       3770, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 16952, 3, 3143,
                                                                       3836, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 17150, 3, 3198,
                                                                       3902, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 17348, 3, 3308,
                                                                       3968, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 17546, 3, 3363,
                                                                       4034, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 17744, 3, 3418,
                                                                       4100, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 17942, 3, 3473,
                                                                       4166, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 18140, 3, 3528,
                                                                       4232, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 18338, 3, 3638,
                                                                       4298, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 18572, 3, 3704,
                                                                       4376, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 18806, 3, 3770,
                                                                       4454, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 19040, 3, 3836,
                                                                       4532, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 19274, 3, 3968,
                                                                       4610, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 19508, 3, 4034,
                                                                       4688, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 19742, 3, 4100,
                                                                       4766, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 19976, 3, 4166,
                                                                       4844, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 20210, 3, 4298,
                                                                       4922, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 20483, 3, 4376,
                                                                       5013, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 20756, 3, 4454,
                                                                       5104, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 21029, 3, 4610,
                                                                       5195, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 21302, 3, 4688,
                                                                       5286, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 21575, 3, 4766,
                                                                       5377, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21848, 3, 7, 8,
                                                                       5474, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21854, 3, 8, 9,
                                                                       5477, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21860, 3, 9, 10,
                                                                       5480, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21866, 3, 10, 11,
                                                                       5483, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21872, 3, 11, 12,
                                                                       5486, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21878, 3, 12, 13,
                                                                       5489, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21884, 3, 13, 14,
                                                                       5492, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21890, 3, 14, 15,
                                                                       5495, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21896, 3, 15, 16,
                                                                       5498, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21902, 3, 16, 17,
                                                                       5501, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21908, 3, 17, 18,
                                                                       5504, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21914, 3, 18, 19,
                                                                       5507, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21920, 3, 19, 20,
                                                                       5510, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21926, 3, 23, 24,
                                                                       5519, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21932, 3, 24, 25,
                                                                       5522, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21938, 3, 25, 26,
                                                                       5525, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21944, 3, 26, 27,
                                                                       5528, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21950, 3, 27, 28,
                                                                       5531, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21956, 3, 28, 29,
                                                                       5534, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21962, 3, 29, 30,
                                                                       5537, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21968, 3, 30, 31,
                                                                       5540, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21974, 3, 31, 32,
                                                                       5543, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21980, 3, 32, 33,
                                                                       5546, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21986, 3, 33, 34,
                                                                       5549, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21992, 3, 34, 35,
                                                                       5552, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21998, 3, 35, 36,
                                                                       5555, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22004, 0, 3,
                                                                       21848, 5474, 21854, 5576,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22022, 0, 3,
                                                                       21854, 5477, 21860, 5585,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22040, 0, 3,
                                                                       21860, 5480, 21866, 5594,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22058, 0, 3,
                                                                       21866, 5483, 21872, 5603,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22076, 0, 3,
                                                                       21872, 5486, 21878, 5612,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22094, 0, 3,
                                                                       21878, 5489, 21884, 5621,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22112, 0, 3,
                                                                       21884, 5492, 21890, 5630,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22130, 0, 3,
                                                                       21890, 5495, 21896, 5639,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22148, 0, 3,
                                                                       21896, 5498, 21902, 5648,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22166, 0, 3,
                                                                       21902, 5501, 21908, 5657,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22184, 0, 3,
                                                                       21908, 5504, 21914, 5666,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22202, 0, 3,
                                                                       21914, 5507, 21920, 5675,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22220, 0, 3,
                                                                       21926, 5519, 21932, 5702,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22238, 0, 3,
                                                                       21932, 5522, 21938, 5711,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22256, 0, 3,
                                                                       21938, 5525, 21944, 5720,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22274, 0, 3,
                                                                       21944, 5528, 21950, 5729,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22292, 0, 3,
                                                                       21950, 5531, 21956, 5738,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22310, 0, 3,
                                                                       21956, 5534, 21962, 5747,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22328, 0, 3,
                                                                       21962, 5537, 21968, 5756,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22346, 0, 3,
                                                                       21968, 5540, 21974, 5765,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22364, 0, 3,
                                                                       21974, 5543, 21980, 5774,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22382, 0, 3,
                                                                       21980, 5546, 21986, 5783,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22400, 0, 3,
                                                                       21986, 5549, 21992, 5792,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 22418, 0, 3,
                                                                       21992, 5552, 21998, 5801,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22436, 0, 3,
                                                                       22004, 5576, 22022, 122,
                                                                       128, 5846, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22472, 0, 3,
                                                                       22022, 5585, 22040, 128,
                                                                       134, 5864, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22508, 0, 3,
                                                                       22040, 5594, 22058, 134,
                                                                       140, 5882, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22544, 0, 3,
                                                                       22058, 5603, 22076, 140,
                                                                       146, 5900, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22580, 0, 3,
                                                                       22076, 5612, 22094, 146,
                                                                       152, 5918, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22616, 0, 3,
                                                                       22094, 5621, 22112, 152,
                                                                       158, 5936, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22652, 0, 3,
                                                                       22112, 5630, 22130, 158,
                                                                       164, 5954, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22688, 0, 3,
                                                                       22130, 5639, 22148, 164,
                                                                       170, 5972, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22724, 0, 3,
                                                                       22148, 5648, 22166, 170,
                                                                       176, 5990, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22760, 0, 3,
                                                                       22166, 5657, 22184, 176,
                                                                       182, 6008, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22796, 0, 3,
                                                                       22184, 5666, 22202, 182,
                                                                       188, 6026, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22832, 0, 3,
                                                                       22220, 5702, 22238, 200,
                                                                       206, 6080, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22868, 0, 3,
                                                                       22238, 5711, 22256, 206,
                                                                       212, 6098, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22904, 0, 3,
                                                                       22256, 5720, 22274, 212,
                                                                       218, 6116, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22940, 0, 3,
                                                                       22274, 5729, 22292, 218,
                                                                       224, 6134, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 22976, 0, 3,
                                                                       22292, 5738, 22310, 224,
                                                                       230, 6152, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23012, 0, 3,
                                                                       22310, 5747, 22328, 230,
                                                                       236, 6170, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23048, 0, 3,
                                                                       22328, 5756, 22346, 236,
                                                                       242, 6188, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23084, 0, 3,
                                                                       22346, 5765, 22364, 242,
                                                                       248, 6206, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23120, 0, 3,
                                                                       22364, 5774, 22382, 248,
                                                                       254, 6224, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23156, 0, 3,
                                                                       22382, 5783, 22400, 254,
                                                                       260, 6242, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 23192, 0, 3,
                                                                       22400, 5792, 22418, 260,
                                                                       266, 6260, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23228, 0, 3,
                                                                       22436, 5846, 22472, 278,
                                                                       288, 6338, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23288, 0, 3,
                                                                       22472, 5864, 22508, 288,
                                                                       298, 6368, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23348, 0, 3,
                                                                       22508, 5882, 22544, 298,
                                                                       308, 6398, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23408, 0, 3,
                                                                       22544, 5900, 22580, 308,
                                                                       318, 6428, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23468, 0, 3,
                                                                       22580, 5918, 22616, 318,
                                                                       328, 6458, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23528, 0, 3,
                                                                       22616, 5936, 22652, 328,
                                                                       338, 6488, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23588, 0, 3,
                                                                       22652, 5954, 22688, 338,
                                                                       348, 6518, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23648, 0, 3,
                                                                       22688, 5972, 22724, 348,
                                                                       358, 6548, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23708, 0, 3,
                                                                       22724, 5990, 22760, 358,
                                                                       368, 6578, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23768, 0, 3,
                                                                       22760, 6008, 22796, 368,
                                                                       378, 6608, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23828, 0, 3,
                                                                       22832, 6080, 22868, 398,
                                                                       408, 6698, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23888, 0, 3,
                                                                       22868, 6098, 22904, 408,
                                                                       418, 6728, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 23948, 0, 3,
                                                                       22904, 6116, 22940, 418,
                                                                       428, 6758, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24008, 0, 3,
                                                                       22940, 6134, 22976, 428,
                                                                       438, 6788, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24068, 0, 3,
                                                                       22976, 6152, 23012, 438,
                                                                       448, 6818, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24128, 0, 3,
                                                                       23012, 6170, 23048, 448,
                                                                       458, 6848, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24188, 0, 3,
                                                                       23048, 6188, 23084, 458,
                                                                       468, 6878, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24248, 0, 3,
                                                                       23084, 6206, 23120, 468,
                                                                       478, 6908, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24308, 0, 3,
                                                                       23120, 6224, 23156, 478,
                                                                       488, 6938, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 24368, 0, 3,
                                                                       23156, 6242, 23192, 488,
                                                                       498, 6968, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 24428, 0, 3,
                                                                       23228, 6338, 23288, 518,
                                                                       533, 7088, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 24518, 0, 3,
                                                                       23288, 6368, 23348, 533,
                                                                       548, 7133, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 24608, 0, 3,
                                                                       23348, 6398, 23408, 548,
                                                                       563, 7178, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 24698, 0, 3,
                                                                       23408, 6428, 23468, 563,
                                                                       578, 7223, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 24788, 0, 3,
                                                                       23468, 6458, 23528, 578,
                                                                       593, 7268, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 24878, 0, 3,
                                                                       23528, 6488, 23588, 593,
                                                                       608, 7313, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 24968, 0, 3,
                                                                       23588, 6518, 23648, 608,
                                                                       623, 7358, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25058, 0, 3,
                                                                       23648, 6548, 23708, 623,
                                                                       638, 7403, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25148, 0, 3,
                                                                       23708, 6578, 23768, 638,
                                                                       653, 7448, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25238, 0, 3,
                                                                       23828, 6698, 23888, 683,
                                                                       698, 7583, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25328, 0, 3,
                                                                       23888, 6728, 23948, 698,
                                                                       713, 7628, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25418, 0, 3,
                                                                       23948, 6758, 24008, 713,
                                                                       728, 7673, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25508, 0, 3,
                                                                       24008, 6788, 24068, 728,
                                                                       743, 7718, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25598, 0, 3,
                                                                       24068, 6818, 24128, 743,
                                                                       758, 7763, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25688, 0, 3,
                                                                       24128, 6848, 24188, 758,
                                                                       773, 7808, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25778, 0, 3,
                                                                       24188, 6878, 24248, 773,
                                                                       788, 7853, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25868, 0, 3,
                                                                       24248, 6908, 24308, 788,
                                                                       803, 7898, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 25958, 0, 3,
                                                                       24308, 6938, 24368, 803,
                                                                       818, 7943, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 26048, 0, 3,
                                                                       24428, 7088, 24518, 848,
                                                                       869, 8114, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 26174, 0, 3,
                                                                       24518, 7133, 24608, 869,
                                                                       890, 8177, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 26300, 0, 3,
                                                                       24608, 7178, 24698, 890,
                                                                       911, 8240, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 26426, 0, 3,
                                                                       24698, 7223, 24788, 911,
                                                                       932, 8303, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 26552, 0, 3,
                                                                       24788, 7268, 24878, 932,
                                                                       953, 8366, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 26678, 0, 3,
                                                                       24878, 7313, 24968, 953,
                                                                       974, 8429, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 26804, 0, 3,
                                                                       24968, 7358, 25058, 974,
                                                                       995, 8492, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 26930, 0, 3,
                                                                       25058, 7403, 25148, 995,
                                                                       1016, 8555, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27056, 0, 3,
                                                                       25238, 7583, 25328, 1058,
                                                                       1079, 8744, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27182, 0, 3,
                                                                       25328, 7628, 25418, 1079,
                                                                       1100, 8807, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27308, 0, 3,
                                                                       25418, 7673, 25508, 1100,
                                                                       1121, 8870, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27434, 0, 3,
                                                                       25508, 7718, 25598, 1121,
                                                                       1142, 8933, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27560, 0, 3,
                                                                       25598, 7763, 25688, 1142,
                                                                       1163, 8996, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27686, 0, 3,
                                                                       25688, 7808, 25778, 1163,
                                                                       1184, 9059, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27812, 0, 3,
                                                                       25778, 7853, 25868, 1184,
                                                                       1205, 9122, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 27938, 0, 3,
                                                                       25868, 7898, 25958, 1205,
                                                                       1226, 9185, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 28064, 0, 3,
                                                                       26048, 8114, 26174, 1268,
                                                                       1296, 9416, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 28232, 0, 3,
                                                                       26174, 8177, 26300, 1296,
                                                                       1324, 9500, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 28400, 0, 3,
                                                                       26300, 8240, 26426, 1324,
                                                                       1352, 9584, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 28568, 0, 3,
                                                                       26426, 8303, 26552, 1352,
                                                                       1380, 9668, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 28736, 0, 3,
                                                                       26552, 8366, 26678, 1380,
                                                                       1408, 9752, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 28904, 0, 3,
                                                                       26678, 8429, 26804, 1408,
                                                                       1436, 9836, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 29072, 0, 3,
                                                                       26804, 8492, 26930, 1436,
                                                                       1464, 9920, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 29240, 0, 3,
                                                                       27056, 8744, 27182, 1520,
                                                                       1548, 10172, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 29408, 0, 3,
                                                                       27182, 8807, 27308, 1548,
                                                                       1576, 10256, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 29576, 0, 3,
                                                                       27308, 8870, 27434, 1576,
                                                                       1604, 10340, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 29744, 0, 3,
                                                                       27434, 8933, 27560, 1604,
                                                                       1632, 10424, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 29912, 0, 3,
                                                                       27560, 8996, 27686, 1632,
                                                                       1660, 10508, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 30080, 0, 3,
                                                                       27686, 9059, 27812, 1660,
                                                                       1688, 10592, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 30248, 0, 3,
                                                                       27812, 9122, 27938, 1688,
                                                                       1716, 10676, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 30416, 0, 3,
                                                                       28064, 9416, 28232, 1772,
                                                                       1808, 10976, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 30632, 0, 3,
                                                                       28232, 9500, 28400, 1808,
                                                                       1844, 11084, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 30848, 0, 3,
                                                                       28400, 9584, 28568, 1844,
                                                                       1880, 11192, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 31064, 0, 3,
                                                                       28568, 9668, 28736, 1880,
                                                                       1916, 11300, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 31280, 0, 3,
                                                                       28736, 9752, 28904, 1916,
                                                                       1952, 11408, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 31496, 0, 3,
                                                                       28904, 9836, 29072, 1952,
                                                                       1988, 11516, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 31712, 0, 3,
                                                                       29240, 10172, 29408, 2060,
                                                                       2096, 11840, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 31928, 0, 3,
                                                                       29408, 10256, 29576, 2096,
                                                                       2132, 11948, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 32144, 0, 3,
                                                                       29576, 10340, 29744, 2132,
                                                                       2168, 12056, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 32360, 0, 3,
                                                                       29744, 10424, 29912, 2168,
                                                                       2204, 12164, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 32576, 0, 3,
                                                                       29912, 10508, 30080, 2204,
                                                                       2240, 12272, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 32792, 0, 3,
                                                                       30080, 10592, 30248, 2240,
                                                                       2276, 12380, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 33008, 0, 3,
                                                                       30416, 10976, 30632, 2348,
                                                                       2393, 12758, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 33278, 0, 3,
                                                                       30632, 11084, 30848, 2393,
                                                                       2438, 12893, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 33548, 0, 3,
                                                                       30848, 11192, 31064, 2438,
                                                                       2483, 13028, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 33818, 0, 3,
                                                                       31064, 11300, 31280, 2483,
                                                                       2528, 13163, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 34088, 0, 3,
                                                                       31280, 11408, 31496, 2528,
                                                                       2573, 13298, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 34358, 0, 3,
                                                                       31712, 11840, 31928, 2663,
                                                                       2708, 13703, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 34628, 0, 3,
                                                                       31928, 11948, 32144, 2708,
                                                                       2753, 13838, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 34898, 0, 3,
                                                                       32144, 12056, 32360, 2753,
                                                                       2798, 13973, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 35168, 0, 3,
                                                                       32360, 12164, 32576, 2798,
                                                                       2843, 14108, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 35438, 0, 3,
                                                                       32576, 12272, 32792, 2843,
                                                                       2888, 14243, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 35708, 0, 3,
                                                                       33008, 12758, 33278, 2978,
                                                                       3033, 14708, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 36038, 0, 3,
                                                                       33278, 12893, 33548, 3033,
                                                                       3088, 14873, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 36368, 0, 3,
                                                                       33548, 13028, 33818, 3088,
                                                                       3143, 15038, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 36698, 0, 3,
                                                                       33818, 13163, 34088, 3143,
                                                                       3198, 15203, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 37028, 0, 3,
                                                                       34358, 13703, 34628, 3308,
                                                                       3363, 15698, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 37358, 0, 3,
                                                                       34628, 13838, 34898, 3363,
                                                                       3418, 15863, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 37688, 0, 3,
                                                                       34898, 13973, 35168, 3418,
                                                                       3473, 16028, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 38018, 0, 3,
                                                                       35168, 14108, 35438, 3473,
                                                                       3528, 16193, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 38348, 0, 3,
                                                                       35708, 14708, 36038, 3638,
                                                                       3704, 16754, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 38744, 0, 3,
                                                                       36038, 14873, 36368, 3704,
                                                                       3770, 16952, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 39140, 0, 3,
                                                                       36368, 15038, 36698, 3770,
                                                                       3836, 17150, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 39536, 0, 3,
                                                                       37028, 15698, 37358, 3968,
                                                                       4034, 17744, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 39932, 0, 3,
                                                                       37358, 15863, 37688, 4034,
                                                                       4100, 17942, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 40328, 0, 3,
                                                                       37688, 16028, 38018, 4100,
                                                                       4166, 18140, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 40724, 0, 3,
                                                                       38348, 16754, 38744, 4298,
                                                                       4376, 18806, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 41192, 0, 3,
                                                                       38744, 16952, 39140, 4376,
                                                                       4454, 19040, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 41660, 0, 3,
                                                                       39536, 17744, 39932, 4610,
                                                                       4688, 19742, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 42128, 0, 3,
                                                                       39932, 17942, 40328, 4688,
                                                                       4766, 19976, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 42596, 0, 3,
                                                                       40724, 18806, 41192, 4922,
                                                                       5013, 20756, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 43142, 0, 3,
                                                                       41660, 19742, 42128, 5195,
                                                                       5286, 21575, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43688, 3, 5468,
                                                                       5471, 21848, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43698, 3, 5471,
                                                                       5474, 21854, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43708, 3, 5474,
                                                                       5477, 21860, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43718, 3, 5477,
                                                                       5480, 21866, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43728, 3, 5480,
                                                                       5483, 21872, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43738, 3, 5483,
                                                                       5486, 21878, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43748, 3, 5486,
                                                                       5489, 21884, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43758, 3, 5489,
                                                                       5492, 21890, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43768, 3, 5492,
                                                                       5495, 21896, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43778, 3, 5495,
                                                                       5498, 21902, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43788, 3, 5498,
                                                                       5501, 21908, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43798, 3, 5501,
                                                                       5504, 21914, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43808, 3, 5504,
                                                                       5507, 21920, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43818, 3, 5513,
                                                                       5516, 21926, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43828, 3, 5516,
                                                                       5519, 21932, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43838, 3, 5519,
                                                                       5522, 21938, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43848, 3, 5522,
                                                                       5525, 21944, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43858, 3, 5525,
                                                                       5528, 21950, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43868, 3, 5528,
                                                                       5531, 21956, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43878, 3, 5531,
                                                                       5534, 21962, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43888, 3, 5534,
                                                                       5537, 21968, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43898, 3, 5537,
                                                                       5540, 21974, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43908, 3, 5540,
                                                                       5543, 21980, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43918, 3, 5543,
                                                                       5546, 21986, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43928, 3, 5546,
                                                                       5549, 21992, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 43938, 3, 5549,
                                                                       5552, 21998, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 43948, 0, 3,
                                                                       43688, 21848, 43698, 5558,
                                                                       5567, 22004, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 43978, 0, 3,
                                                                       43698, 21854, 43708, 5567,
                                                                       5576, 22022, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44008, 0, 3,
                                                                       43708, 21860, 43718, 5576,
                                                                       5585, 22040, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44038, 0, 3,
                                                                       43718, 21866, 43728, 5585,
                                                                       5594, 22058, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44068, 0, 3,
                                                                       43728, 21872, 43738, 5594,
                                                                       5603, 22076, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44098, 0, 3,
                                                                       43738, 21878, 43748, 5603,
                                                                       5612, 22094, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44128, 0, 3,
                                                                       43748, 21884, 43758, 5612,
                                                                       5621, 22112, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44158, 0, 3,
                                                                       43758, 21890, 43768, 5621,
                                                                       5630, 22130, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44188, 0, 3,
                                                                       43768, 21896, 43778, 5630,
                                                                       5639, 22148, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44218, 0, 3,
                                                                       43778, 21902, 43788, 5639,
                                                                       5648, 22166, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44248, 0, 3,
                                                                       43788, 21908, 43798, 5648,
                                                                       5657, 22184, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44278, 0, 3,
                                                                       43798, 21914, 43808, 5657,
                                                                       5666, 22202, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44308, 0, 3,
                                                                       43818, 21926, 43828, 5684,
                                                                       5693, 22220, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44338, 0, 3,
                                                                       43828, 21932, 43838, 5693,
                                                                       5702, 22238, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44368, 0, 3,
                                                                       43838, 21938, 43848, 5702,
                                                                       5711, 22256, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44398, 0, 3,
                                                                       43848, 21944, 43858, 5711,
                                                                       5720, 22274, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44428, 0, 3,
                                                                       43858, 21950, 43868, 5720,
                                                                       5729, 22292, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44458, 0, 3,
                                                                       43868, 21956, 43878, 5729,
                                                                       5738, 22310, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44488, 0, 3,
                                                                       43878, 21962, 43888, 5738,
                                                                       5747, 22328, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44518, 0, 3,
                                                                       43888, 21968, 43898, 5747,
                                                                       5756, 22346, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44548, 0, 3,
                                                                       43898, 21974, 43908, 5756,
                                                                       5765, 22364, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44578, 0, 3,
                                                                       43908, 21980, 43918, 5765,
                                                                       5774, 22382, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44608, 0, 3,
                                                                       43918, 21986, 43928, 5774,
                                                                       5783, 22400, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 44638, 0, 3,
                                                                       43928, 21992, 43938, 5783,
                                                                       5792, 22418, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 44668, 0, 3,
                                                                       43948, 22004, 43978, 5810,
                                                                       5828, 22436, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 44728, 0, 3,
                                                                       43978, 22022, 44008, 5828,
                                                                       5846, 22472, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 44788, 0, 3,
                                                                       44008, 22040, 44038, 5846,
                                                                       5864, 22508, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 44848, 0, 3,
                                                                       44038, 22058, 44068, 5864,
                                                                       5882, 22544, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 44908, 0, 3,
                                                                       44068, 22076, 44098, 5882,
                                                                       5900, 22580, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 44968, 0, 3,
                                                                       44098, 22094, 44128, 5900,
                                                                       5918, 22616, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45028, 0, 3,
                                                                       44128, 22112, 44158, 5918,
                                                                       5936, 22652, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45088, 0, 3,
                                                                       44158, 22130, 44188, 5936,
                                                                       5954, 22688, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45148, 0, 3,
                                                                       44188, 22148, 44218, 5954,
                                                                       5972, 22724, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45208, 0, 3,
                                                                       44218, 22166, 44248, 5972,
                                                                       5990, 22760, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45268, 0, 3,
                                                                       44248, 22184, 44278, 5990,
                                                                       6008, 22796, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45328, 0, 3,
                                                                       44308, 22220, 44338, 6044,
                                                                       6062, 22832, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45388, 0, 3,
                                                                       44338, 22238, 44368, 6062,
                                                                       6080, 22868, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45448, 0, 3,
                                                                       44368, 22256, 44398, 6080,
                                                                       6098, 22904, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45508, 0, 3,
                                                                       44398, 22274, 44428, 6098,
                                                                       6116, 22940, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45568, 0, 3,
                                                                       44428, 22292, 44458, 6116,
                                                                       6134, 22976, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45628, 0, 3,
                                                                       44458, 22310, 44488, 6134,
                                                                       6152, 23012, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45688, 0, 3,
                                                                       44488, 22328, 44518, 6152,
                                                                       6170, 23048, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45748, 0, 3,
                                                                       44518, 22346, 44548, 6170,
                                                                       6188, 23084, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45808, 0, 3,
                                                                       44548, 22364, 44578, 6188,
                                                                       6206, 23120, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45868, 0, 3,
                                                                       44578, 22382, 44608, 6206,
                                                                       6224, 23156, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 45928, 0, 3,
                                                                       44608, 22400, 44638, 6224,
                                                                       6242, 23192, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 45988, 0, 3,
                                                                       44668, 22436, 44728, 6278,
                                                                       6308, 23228, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46088, 0, 3,
                                                                       44728, 22472, 44788, 6308,
                                                                       6338, 23288, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46188, 0, 3,
                                                                       44788, 22508, 44848, 6338,
                                                                       6368, 23348, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46288, 0, 3,
                                                                       44848, 22544, 44908, 6368,
                                                                       6398, 23408, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46388, 0, 3,
                                                                       44908, 22580, 44968, 6398,
                                                                       6428, 23468, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46488, 0, 3,
                                                                       44968, 22616, 45028, 6428,
                                                                       6458, 23528, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46588, 0, 3,
                                                                       45028, 22652, 45088, 6458,
                                                                       6488, 23588, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46688, 0, 3,
                                                                       45088, 22688, 45148, 6488,
                                                                       6518, 23648, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46788, 0, 3,
                                                                       45148, 22724, 45208, 6518,
                                                                       6548, 23708, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46888, 0, 3,
                                                                       45208, 22760, 45268, 6548,
                                                                       6578, 23768, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 46988, 0, 3,
                                                                       45328, 22832, 45388, 6638,
                                                                       6668, 23828, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47088, 0, 3,
                                                                       45388, 22868, 45448, 6668,
                                                                       6698, 23888, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47188, 0, 3,
                                                                       45448, 22904, 45508, 6698,
                                                                       6728, 23948, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47288, 0, 3,
                                                                       45508, 22940, 45568, 6728,
                                                                       6758, 24008, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47388, 0, 3,
                                                                       45568, 22976, 45628, 6758,
                                                                       6788, 24068, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47488, 0, 3,
                                                                       45628, 23012, 45688, 6788,
                                                                       6818, 24128, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47588, 0, 3,
                                                                       45688, 23048, 45748, 6818,
                                                                       6848, 24188, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47688, 0, 3,
                                                                       45748, 23084, 45808, 6848,
                                                                       6878, 24248, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47788, 0, 3,
                                                                       45808, 23120, 45868, 6878,
                                                                       6908, 24308, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47888, 0, 3,
                                                                       45868, 23156, 45928, 6908,
                                                                       6938, 24368, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 47988, 0, 3,
                                                                       45988, 23228, 46088, 6998,
                                                                       7043, 24428, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48138, 0, 3,
                                                                       46088, 23288, 46188, 7043,
                                                                       7088, 24518, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48288, 0, 3,
                                                                       46188, 23348, 46288, 7088,
                                                                       7133, 24608, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48438, 0, 3,
                                                                       46288, 23408, 46388, 7133,
                                                                       7178, 24698, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48588, 0, 3,
                                                                       46388, 23468, 46488, 7178,
                                                                       7223, 24788, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48738, 0, 3,
                                                                       46488, 23528, 46588, 7223,
                                                                       7268, 24878, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 48888, 0, 3,
                                                                       46588, 23588, 46688, 7268,
                                                                       7313, 24968, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49038, 0, 3,
                                                                       46688, 23648, 46788, 7313,
                                                                       7358, 25058, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49188, 0, 3,
                                                                       46788, 23708, 46888, 7358,
                                                                       7403, 25148, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49338, 0, 3,
                                                                       46988, 23828, 47088, 7493,
                                                                       7538, 25238, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49488, 0, 3,
                                                                       47088, 23888, 47188, 7538,
                                                                       7583, 25328, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49638, 0, 3,
                                                                       47188, 23948, 47288, 7583,
                                                                       7628, 25418, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49788, 0, 3,
                                                                       47288, 24008, 47388, 7628,
                                                                       7673, 25508, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49938, 0, 3,
                                                                       47388, 24068, 47488, 7673,
                                                                       7718, 25598, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50088, 0, 3,
                                                                       47488, 24128, 47588, 7718,
                                                                       7763, 25688, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50238, 0, 3,
                                                                       47588, 24188, 47688, 7763,
                                                                       7808, 25778, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50388, 0, 3,
                                                                       47688, 24248, 47788, 7808,
                                                                       7853, 25868, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50538, 0, 3,
                                                                       47788, 24308, 47888, 7853,
                                                                       7898, 25958, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 50688, 0, 3,
                                                                       47988, 24428, 48138, 7988,
                                                                       8051, 26048, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 50898, 0, 3,
                                                                       48138, 24518, 48288, 8051,
                                                                       8114, 26174, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51108, 0, 3,
                                                                       48288, 24608, 48438, 8114,
                                                                       8177, 26300, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51318, 0, 3,
                                                                       48438, 24698, 48588, 8177,
                                                                       8240, 26426, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51528, 0, 3,
                                                                       48588, 24788, 48738, 8240,
                                                                       8303, 26552, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51738, 0, 3,
                                                                       48738, 24878, 48888, 8303,
                                                                       8366, 26678, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51948, 0, 3,
                                                                       48888, 24968, 49038, 8366,
                                                                       8429, 26804, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52158, 0, 3,
                                                                       49038, 25058, 49188, 8429,
                                                                       8492, 26930, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52368, 0, 3,
                                                                       49338, 25238, 49488, 8618,
                                                                       8681, 27056, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52578, 0, 3,
                                                                       49488, 25328, 49638, 8681,
                                                                       8744, 27182, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52788, 0, 3,
                                                                       49638, 25418, 49788, 8744,
                                                                       8807, 27308, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52998, 0, 3,
                                                                       49788, 25508, 49938, 8807,
                                                                       8870, 27434, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 53208, 0, 3,
                                                                       49938, 25598, 50088, 8870,
                                                                       8933, 27560, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 53418, 0, 3,
                                                                       50088, 25688, 50238, 8933,
                                                                       8996, 27686, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 53628, 0, 3,
                                                                       50238, 25778, 50388, 8996,
                                                                       9059, 27812, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 53838, 0, 3,
                                                                       50388, 25868, 50538, 9059,
                                                                       9122, 27938, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54048, 0, 3,
                                                                       50688, 26048, 50898, 9248,
                                                                       9332, 28064, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54328, 0, 3,
                                                                       50898, 26174, 51108, 9332,
                                                                       9416, 28232, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54608, 0, 3,
                                                                       51108, 26300, 51318, 9416,
                                                                       9500, 28400, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54888, 0, 3,
                                                                       51318, 26426, 51528, 9500,
                                                                       9584, 28568, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 55168, 0, 3,
                                                                       51528, 26552, 51738, 9584,
                                                                       9668, 28736, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 55448, 0, 3,
                                                                       51738, 26678, 51948, 9668,
                                                                       9752, 28904, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 55728, 0, 3,
                                                                       51948, 26804, 52158, 9752,
                                                                       9836, 29072, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 56008, 0, 3,
                                                                       52368, 27056, 52578,
                                                                       10004, 10088, 29240,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 56288, 0, 3,
                                                                       52578, 27182, 52788,
                                                                       10088, 10172, 29408,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 56568, 0, 3,
                                                                       52788, 27308, 52998,
                                                                       10172, 10256, 29576,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 56848, 0, 3,
                                                                       52998, 27434, 53208,
                                                                       10256, 10340, 29744,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 57128, 0, 3,
                                                                       53208, 27560, 53418,
                                                                       10340, 10424, 29912,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 57408, 0, 3,
                                                                       53418, 27686, 53628,
                                                                       10424, 10508, 30080,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 57688, 0, 3,
                                                                       53628, 27812, 53838,
                                                                       10508, 10592, 30248,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 57968, 0, 3,
                                                                       54048, 28064, 54328,
                                                                       10760, 10868, 30416,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 58328, 0, 3,
                                                                       54328, 28232, 54608,
                                                                       10868, 10976, 30632,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 58688, 0, 3,
                                                                       54608, 28400, 54888,
                                                                       10976, 11084, 30848,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 59048, 0, 3,
                                                                       54888, 28568, 55168,
                                                                       11084, 11192, 31064,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 59408, 0, 3,
                                                                       55168, 28736, 55448,
                                                                       11192, 11300, 31280,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 59768, 0, 3,
                                                                       55448, 28904, 55728,
                                                                       11300, 11408, 31496,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 60128, 0, 3,
                                                                       56008, 29240, 56288,
                                                                       11624, 11732, 31712,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 60488, 0, 3,
                                                                       56288, 29408, 56568,
                                                                       11732, 11840, 31928,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 60848, 0, 3,
                                                                       56568, 29576, 56848,
                                                                       11840, 11948, 32144,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 61208, 0, 3,
                                                                       56848, 29744, 57128,
                                                                       11948, 12056, 32360,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 61568, 0, 3,
                                                                       57128, 29912, 57408,
                                                                       12056, 12164, 32576,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 61928, 0, 3,
                                                                       57408, 30080, 57688,
                                                                       12164, 12272, 32792,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 62288, 0, 3,
                                                                       57968, 30416, 58328,
                                                                       12488, 12623, 33008,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 62738, 0, 3,
                                                                       58328, 30632, 58688,
                                                                       12623, 12758, 33278,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 63188, 0, 3,
                                                                       58688, 30848, 59048,
                                                                       12758, 12893, 33548,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 63638, 0, 3,
                                                                       59048, 31064, 59408,
                                                                       12893, 13028, 33818,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 64088, 0, 3,
                                                                       59408, 31280, 59768,
                                                                       13028, 13163, 34088,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 64538, 0, 3,
                                                                       60128, 31712, 60488,
                                                                       13433, 13568, 34358,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 64988, 0, 3,
                                                                       60488, 31928, 60848,
                                                                       13568, 13703, 34628,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 65438, 0, 3,
                                                                       60848, 32144, 61208,
                                                                       13703, 13838, 34898,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 65888, 0, 3,
                                                                       61208, 32360, 61568,
                                                                       13838, 13973, 35168,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 66338, 0, 3,
                                                                       61568, 32576, 61928,
                                                                       13973, 14108, 35438,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 66788, 0, 3,
                                                                       62288, 33008, 62738,
                                                                       14378, 14543, 35708,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 67338, 0, 3,
                                                                       62738, 33278, 63188,
                                                                       14543, 14708, 36038,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 67888, 0, 3,
                                                                       63188, 33548, 63638,
                                                                       14708, 14873, 36368,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 68438, 0, 3,
                                                                       63638, 33818, 64088,
                                                                       14873, 15038, 36698,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 68988, 0, 3,
                                                                       64538, 34358, 64988,
                                                                       15368, 15533, 37028,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 69538, 0, 3,
                                                                       64988, 34628, 65438,
                                                                       15533, 15698, 37358,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 70088, 0, 3,
                                                                       65438, 34898, 65888,
                                                                       15698, 15863, 37688,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 70638, 0, 3,
                                                                       65888, 35168, 66338,
                                                                       15863, 16028, 38018,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 71188, 0, 3,
                                                                       66788, 35708, 67338,
                                                                       16358, 16556, 38348,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 71848, 0, 3,
                                                                       67338, 36038, 67888,
                                                                       16556, 16754, 38744,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 72508, 0, 3,
                                                                       67888, 36368, 68438,
                                                                       16754, 16952, 39140,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 73168, 0, 3,
                                                                       68988, 37028, 69538,
                                                                       17348, 17546, 39536,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 73828, 0, 3,
                                                                       69538, 37358, 70088,
                                                                       17546, 17744, 39932,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 74488, 0, 3,
                                                                       70088, 37688, 70638,
                                                                       17744, 17942, 40328,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 75148, 0, 3,
                                                                       71188, 38348, 71848,
                                                                       18338, 18572, 40724,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 75928, 0, 3,
                                                                       71848, 38744, 72508,
                                                                       18572, 18806, 41192,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 76708, 0, 3,
                                                                       73168, 39536, 73828,
                                                                       19274, 19508, 41660,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 77488, 0, 3,
                                                                       73828, 39932, 74488,
                                                                       19508, 19742, 42128,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 78268, 0, 3,
                                                                       75148, 40724, 75928,
                                                                       20210, 20483, 42596,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 79178, 0, 3,
                                                                       76708, 41660, 77488,
                                                                       21029, 21302, 43142,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 80088, 54048, 280, ncols);

                    simdfunc::contract_primitives(buffer, 80564, 56008, 280, ncols);

                    simdfunc::contract_primitives(buffer, 81040, 57968, 360, ncols);

                    simdfunc::contract_primitives(buffer, 81652, 60128, 360, ncols);

                    simdfunc::contract_primitives(buffer, 82264, 62288, 450, ncols);

                    simdfunc::contract_primitives(buffer, 83029, 64538, 450, ncols);

                    simdfunc::contract_primitives(buffer, 83794, 66788, 550, ncols);

                    simdfunc::contract_primitives(buffer, 84729, 68988, 550, ncols);

                    simdfunc::contract_primitives(buffer, 85664, 71188, 660, ncols);

                    simdfunc::contract_primitives(buffer, 86786, 73168, 660, ncols);

                    simdfunc::contract_primitives(buffer, 87908, 75148, 780, ncols);

                    simdfunc::contract_primitives(buffer, 89234, 76708, 780, ncols);

                    simdfunc::contract_primitives(buffer, 90560, 78268, 910, ncols);

                    simdfunc::contract_primitives(buffer, 92107, 79178, 910, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 80368, 80088, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 80844, 80564, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 81400, 81040, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 82012, 81652, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 82714, 82264, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 83479, 83029, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 84344, 83794, 55, 1, nmax);

        simdtrf::transform_f_inner(buffer, 85279, 84729, 55, 1, nmax);

        simdtrf::transform_f_inner(buffer, 86324, 85664, 66, 1, nmax);

        simdtrf::transform_f_inner(buffer, 87446, 86786, 66, 1, nmax);

        simdtrf::transform_f_inner(buffer, 88688, 87908, 78, 1, nmax);

        simdtrf::transform_f_inner(buffer, 90014, 89234, 78, 1, nmax);

        simdtrf::transform_f_inner(buffer, 91470, 90560, 91, 1, nmax);

        simdtrf::transform_f_inner(buffer, 93017, 92107, 91, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 93654, 80368, 81400, 7, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 94242, 80844, 82012, 7, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 94830, 81400, 82714, 7, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 95586, 82012, 83479, 7, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 96342, 82714, 84344, 7, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 97287, 83479, 85279, 7, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 98232, 84344, 86324, 7, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 99387, 85279, 87446, 7, nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 100542, 86324, 88688, 7,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 101928, 87446, 90014, 7,
                                             nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 103314, 88688, 91470, 7,
                                             nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 104952, 90014, 93017, 7,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 106590, 93654, 94830, 7,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 107766, 94242, 95586, 7,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 108942, 94830, 96342, 7,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 110454, 95586, 97287, 7,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 111966, 96342, 98232, 7,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 113856, 97287, 99387, 7,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 115746, 98232, 100542, 7,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 118056, 99387, 101928, 7,
                                             nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 120366, 100542, 103314, 7,
                                             nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 123138, 101928, 104952, 7,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 125910, 106590, 108942, 7,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 127870, 107766, 110454, 7,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 129830, 108942, 111966, 7,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 132350, 110454, 113856, 7,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 134870, 111966, 115746, 7,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 138020, 113856, 118056, 7,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 141170, 115746, 120366, 7,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 145020, 118056, 123138, 7,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 148870, 125910, 129830, 7,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 151810, 127870, 132350, 7,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 154750, 129830, 134870, 7,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 158530, 132350, 138020, 7,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 162310, 134870, 141170, 7,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 167035, 138020, 145020, 7,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 171760, 148870, 154750, 7,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 175876, 151810, 158530, 7,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 179992, 154750, 162310, 7,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 185284, 158530, 167035, 7,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 190576, 171760, 179992, 7,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 196064, 175876, 185284, 7,
                                             nmax);

        simdtrf::transform_i_inner(buffer, 201552, 196064, 28, 7, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 201552, 91, nmax);

        simdtrf::transform_i_inner(buffer, 201552, 190576, 28, 7, nmax);

        simdtrf::transform_i_outer(values + 1183 * nvalues + n * npairs, nvalues, buffer, 201552,
                                   91, nmax);
    }

    for (size_t m = 0; m < 2366; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
