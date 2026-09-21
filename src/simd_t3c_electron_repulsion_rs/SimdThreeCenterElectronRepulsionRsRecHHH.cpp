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


#include "SimdThreeCenterElectronRepulsionRsRecHHH.hpp"

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
#include "SimdTransferHD.hpp"
#include "SimdTransferHF.hpp"
#include "SimdTransferHG.hpp"
#include "SimdTransferHH.hpp"
#include "SimdTransferHP.hpp"
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
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_hhh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_hhh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 232123, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2662 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 232123, 139576, 15338, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4298, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4301, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4304, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4307, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4310, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4313, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4316, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4319, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4322, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4325, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4328, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4331, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4334, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4337, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4340, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4343, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4346, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4349, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4352, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4355, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4358, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4361, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4364, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4367, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4370, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4373, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4376, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4379, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4382, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4385, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4388, 3, 7, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4397, 3, 8, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4406, 3, 9, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4415, 3, 10, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4424, 3, 11, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4433, 3, 12, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4442, 3, 13, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4451, 3, 14, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4460, 3, 15, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4469, 3, 16, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4478, 3, 17, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4487, 3, 18, 71,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4496, 3, 19, 74,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4505, 3, 20, 77,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4514, 3, 23, 80,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4523, 3, 24, 83,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4532, 3, 25, 86,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4541, 3, 26, 89,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4550, 3, 27, 92,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4559, 3, 28, 95,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4568, 3, 29, 98,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4577, 3, 30, 101,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4586, 3, 31, 104,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4595, 3, 32, 107,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4604, 3, 33, 110,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4613, 3, 34, 113,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4622, 3, 35, 116,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4631, 3, 36, 119,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4640, 3, 38, 122,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4658, 3, 41, 128,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4676, 3, 44, 134,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4694, 3, 47, 140,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4712, 3, 50, 146,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4730, 3, 53, 152,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4748, 3, 56, 158,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4766, 3, 59, 164,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4784, 3, 62, 170,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4802, 3, 65, 176,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4820, 3, 68, 182,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4838, 3, 71, 188,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4856, 3, 74, 194,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4874, 3, 80, 200,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4892, 3, 83, 206,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4910, 3, 86, 212,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4928, 3, 89, 218,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4946, 3, 92, 224,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4964, 3, 95, 230,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4982, 3, 98, 236,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5000, 3, 101, 242,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5018, 3, 104, 248,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5036, 3, 107, 254,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5054, 3, 110, 260,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5072, 3, 113, 266,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5090, 3, 116, 272,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5108, 3, 122, 278,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5138, 3, 128, 288,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5168, 3, 134, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5198, 3, 140, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5228, 3, 146, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5258, 3, 152, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5288, 3, 158, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5318, 3, 164, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5348, 3, 170, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5378, 3, 176, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5408, 3, 182, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5438, 3, 188, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5468, 3, 200, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5498, 3, 206, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5528, 3, 212, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5558, 3, 218, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5588, 3, 224, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5618, 3, 230, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5648, 3, 236, 458,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5678, 3, 242, 468,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5708, 3, 248, 478,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5738, 3, 254, 488,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5768, 3, 260, 498,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5798, 3, 266, 508,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5828, 3, 278, 518,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5873, 3, 288, 533,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5918, 3, 298, 548,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5963, 3, 308, 563,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6008, 3, 318, 578,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6053, 3, 328, 593,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6098, 3, 338, 608,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6143, 3, 348, 623,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6188, 3, 358, 638,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6233, 3, 368, 653,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6278, 3, 378, 668,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6323, 3, 398, 683,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6368, 3, 408, 698,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6413, 3, 418, 713,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6458, 3, 428, 728,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6503, 3, 438, 743,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6548, 3, 448, 758,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6593, 3, 458, 773,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6638, 3, 468, 788,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6683, 3, 478, 803,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6728, 3, 488, 818,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6773, 3, 498, 833,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6818, 3, 518, 848,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6881, 3, 533, 869,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6944, 3, 548, 890,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7007, 3, 563, 911,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7070, 3, 578, 932,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7133, 3, 593, 953,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7196, 3, 608, 974,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7259, 3, 623, 995,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7322, 3, 638,
                                                                       1016, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7385, 3, 653,
                                                                       1037, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7448, 3, 683,
                                                                       1058, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7511, 3, 698,
                                                                       1079, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7574, 3, 713,
                                                                       1100, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7637, 3, 728,
                                                                       1121, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7700, 3, 743,
                                                                       1142, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7763, 3, 758,
                                                                       1163, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7826, 3, 773,
                                                                       1184, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7889, 3, 788,
                                                                       1205, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7952, 3, 803,
                                                                       1226, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8015, 3, 818,
                                                                       1247, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8078, 3, 848,
                                                                       1268, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8162, 3, 869,
                                                                       1296, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8246, 3, 890,
                                                                       1324, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8330, 3, 911,
                                                                       1352, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8414, 3, 932,
                                                                       1380, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8498, 3, 953,
                                                                       1408, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8582, 3, 974,
                                                                       1436, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8666, 3, 995,
                                                                       1464, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8750, 3, 1016,
                                                                       1492, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8834, 3, 1058,
                                                                       1520, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8918, 3, 1079,
                                                                       1548, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9002, 3, 1100,
                                                                       1576, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9086, 3, 1121,
                                                                       1604, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9170, 3, 1142,
                                                                       1632, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9254, 3, 1163,
                                                                       1660, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9338, 3, 1184,
                                                                       1688, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9422, 3, 1205,
                                                                       1716, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9506, 3, 1226,
                                                                       1744, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9590, 3, 1268,
                                                                       1772, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9698, 3, 1296,
                                                                       1808, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9806, 3, 1324,
                                                                       1844, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9914, 3, 1352,
                                                                       1880, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10022, 3, 1380,
                                                                       1916, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10130, 3, 1408,
                                                                       1952, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10238, 3, 1436,
                                                                       1988, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10346, 3, 1464,
                                                                       2024, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10454, 3, 1520,
                                                                       2060, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10562, 3, 1548,
                                                                       2096, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10670, 3, 1576,
                                                                       2132, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10778, 3, 1604,
                                                                       2168, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10886, 3, 1632,
                                                                       2204, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10994, 3, 1660,
                                                                       2240, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11102, 3, 1688,
                                                                       2276, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 11210, 3, 1716,
                                                                       2312, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11318, 3, 1772,
                                                                       2348, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11453, 3, 1808,
                                                                       2393, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11588, 3, 1844,
                                                                       2438, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11723, 3, 1880,
                                                                       2483, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11858, 3, 1916,
                                                                       2528, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11993, 3, 1952,
                                                                       2573, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12128, 3, 1988,
                                                                       2618, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12263, 3, 2060,
                                                                       2663, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12398, 3, 2096,
                                                                       2708, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12533, 3, 2132,
                                                                       2753, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12668, 3, 2168,
                                                                       2798, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12803, 3, 2204,
                                                                       2843, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12938, 3, 2240,
                                                                       2888, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 13073, 3, 2276,
                                                                       2933, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13208, 3, 2348,
                                                                       2978, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13373, 3, 2393,
                                                                       3033, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13538, 3, 2438,
                                                                       3088, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13703, 3, 2483,
                                                                       3143, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13868, 3, 2528,
                                                                       3198, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14033, 3, 2573,
                                                                       3253, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14198, 3, 2663,
                                                                       3308, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14363, 3, 2708,
                                                                       3363, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14528, 3, 2753,
                                                                       3418, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14693, 3, 2798,
                                                                       3473, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 14858, 3, 2843,
                                                                       3528, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 15023, 3, 2888,
                                                                       3583, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 15188, 3, 2978,
                                                                       3638, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 15386, 3, 3033,
                                                                       3704, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 15584, 3, 3088,
                                                                       3770, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 15782, 3, 3143,
                                                                       3836, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 15980, 3, 3198,
                                                                       3902, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 16178, 3, 3308,
                                                                       3968, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 16376, 3, 3363,
                                                                       4034, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 16574, 3, 3418,
                                                                       4100, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 16772, 3, 3473,
                                                                       4166, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 16970, 3, 3528,
                                                                       4232, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17168, 3, 7, 8,
                                                                       4304, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17174, 3, 8, 9,
                                                                       4307, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17180, 3, 9, 10,
                                                                       4310, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17186, 3, 10, 11,
                                                                       4313, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17192, 3, 11, 12,
                                                                       4316, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17198, 3, 12, 13,
                                                                       4319, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17204, 3, 13, 14,
                                                                       4322, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17210, 3, 14, 15,
                                                                       4325, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17216, 3, 15, 16,
                                                                       4328, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17222, 3, 16, 17,
                                                                       4331, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17228, 3, 17, 18,
                                                                       4334, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17234, 3, 18, 19,
                                                                       4337, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17240, 3, 19, 20,
                                                                       4340, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17246, 3, 23, 24,
                                                                       4349, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17252, 3, 24, 25,
                                                                       4352, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17258, 3, 25, 26,
                                                                       4355, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17264, 3, 26, 27,
                                                                       4358, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17270, 3, 27, 28,
                                                                       4361, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17276, 3, 28, 29,
                                                                       4364, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17282, 3, 29, 30,
                                                                       4367, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17288, 3, 30, 31,
                                                                       4370, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17294, 3, 31, 32,
                                                                       4373, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17300, 3, 32, 33,
                                                                       4376, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17306, 3, 33, 34,
                                                                       4379, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17312, 3, 34, 35,
                                                                       4382, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 17318, 3, 35, 36,
                                                                       4385, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17324, 0, 3,
                                                                       17168, 4304, 17174, 4406,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17342, 0, 3,
                                                                       17174, 4307, 17180, 4415,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17360, 0, 3,
                                                                       17180, 4310, 17186, 4424,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17378, 0, 3,
                                                                       17186, 4313, 17192, 4433,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17396, 0, 3,
                                                                       17192, 4316, 17198, 4442,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17414, 0, 3,
                                                                       17198, 4319, 17204, 4451,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17432, 0, 3,
                                                                       17204, 4322, 17210, 4460,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17450, 0, 3,
                                                                       17210, 4325, 17216, 4469,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17468, 0, 3,
                                                                       17216, 4328, 17222, 4478,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17486, 0, 3,
                                                                       17222, 4331, 17228, 4487,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17504, 0, 3,
                                                                       17228, 4334, 17234, 4496,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17522, 0, 3,
                                                                       17234, 4337, 17240, 4505,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17540, 0, 3,
                                                                       17246, 4349, 17252, 4532,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17558, 0, 3,
                                                                       17252, 4352, 17258, 4541,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17576, 0, 3,
                                                                       17258, 4355, 17264, 4550,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17594, 0, 3,
                                                                       17264, 4358, 17270, 4559,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17612, 0, 3,
                                                                       17270, 4361, 17276, 4568,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17630, 0, 3,
                                                                       17276, 4364, 17282, 4577,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17648, 0, 3,
                                                                       17282, 4367, 17288, 4586,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17666, 0, 3,
                                                                       17288, 4370, 17294, 4595,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17684, 0, 3,
                                                                       17294, 4373, 17300, 4604,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17702, 0, 3,
                                                                       17300, 4376, 17306, 4613,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17720, 0, 3,
                                                                       17306, 4379, 17312, 4622,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 17738, 0, 3,
                                                                       17312, 4382, 17318, 4631,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17756, 0, 3,
                                                                       17324, 4406, 17342, 122,
                                                                       128, 4676, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17792, 0, 3,
                                                                       17342, 4415, 17360, 128,
                                                                       134, 4694, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17828, 0, 3,
                                                                       17360, 4424, 17378, 134,
                                                                       140, 4712, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17864, 0, 3,
                                                                       17378, 4433, 17396, 140,
                                                                       146, 4730, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17900, 0, 3,
                                                                       17396, 4442, 17414, 146,
                                                                       152, 4748, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17936, 0, 3,
                                                                       17414, 4451, 17432, 152,
                                                                       158, 4766, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17972, 0, 3,
                                                                       17432, 4460, 17450, 158,
                                                                       164, 4784, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18008, 0, 3,
                                                                       17450, 4469, 17468, 164,
                                                                       170, 4802, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18044, 0, 3,
                                                                       17468, 4478, 17486, 170,
                                                                       176, 4820, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18080, 0, 3,
                                                                       17486, 4487, 17504, 176,
                                                                       182, 4838, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18116, 0, 3,
                                                                       17504, 4496, 17522, 182,
                                                                       188, 4856, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18152, 0, 3,
                                                                       17540, 4532, 17558, 200,
                                                                       206, 4910, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18188, 0, 3,
                                                                       17558, 4541, 17576, 206,
                                                                       212, 4928, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18224, 0, 3,
                                                                       17576, 4550, 17594, 212,
                                                                       218, 4946, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18260, 0, 3,
                                                                       17594, 4559, 17612, 218,
                                                                       224, 4964, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18296, 0, 3,
                                                                       17612, 4568, 17630, 224,
                                                                       230, 4982, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18332, 0, 3,
                                                                       17630, 4577, 17648, 230,
                                                                       236, 5000, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18368, 0, 3,
                                                                       17648, 4586, 17666, 236,
                                                                       242, 5018, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18404, 0, 3,
                                                                       17666, 4595, 17684, 242,
                                                                       248, 5036, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18440, 0, 3,
                                                                       17684, 4604, 17702, 248,
                                                                       254, 5054, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18476, 0, 3,
                                                                       17702, 4613, 17720, 254,
                                                                       260, 5072, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18512, 0, 3,
                                                                       17720, 4622, 17738, 260,
                                                                       266, 5090, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18548, 0, 3,
                                                                       17756, 4676, 17792, 278,
                                                                       288, 5168, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18608, 0, 3,
                                                                       17792, 4694, 17828, 288,
                                                                       298, 5198, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18668, 0, 3,
                                                                       17828, 4712, 17864, 298,
                                                                       308, 5228, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18728, 0, 3,
                                                                       17864, 4730, 17900, 308,
                                                                       318, 5258, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18788, 0, 3,
                                                                       17900, 4748, 17936, 318,
                                                                       328, 5288, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18848, 0, 3,
                                                                       17936, 4766, 17972, 328,
                                                                       338, 5318, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18908, 0, 3,
                                                                       17972, 4784, 18008, 338,
                                                                       348, 5348, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18968, 0, 3,
                                                                       18008, 4802, 18044, 348,
                                                                       358, 5378, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19028, 0, 3,
                                                                       18044, 4820, 18080, 358,
                                                                       368, 5408, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19088, 0, 3,
                                                                       18080, 4838, 18116, 368,
                                                                       378, 5438, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19148, 0, 3,
                                                                       18152, 4910, 18188, 398,
                                                                       408, 5528, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19208, 0, 3,
                                                                       18188, 4928, 18224, 408,
                                                                       418, 5558, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19268, 0, 3,
                                                                       18224, 4946, 18260, 418,
                                                                       428, 5588, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19328, 0, 3,
                                                                       18260, 4964, 18296, 428,
                                                                       438, 5618, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19388, 0, 3,
                                                                       18296, 4982, 18332, 438,
                                                                       448, 5648, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19448, 0, 3,
                                                                       18332, 5000, 18368, 448,
                                                                       458, 5678, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19508, 0, 3,
                                                                       18368, 5018, 18404, 458,
                                                                       468, 5708, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19568, 0, 3,
                                                                       18404, 5036, 18440, 468,
                                                                       478, 5738, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19628, 0, 3,
                                                                       18440, 5054, 18476, 478,
                                                                       488, 5768, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19688, 0, 3,
                                                                       18476, 5072, 18512, 488,
                                                                       498, 5798, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19748, 0, 3,
                                                                       18548, 5168, 18608, 518,
                                                                       533, 5918, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19838, 0, 3,
                                                                       18608, 5198, 18668, 533,
                                                                       548, 5963, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19928, 0, 3,
                                                                       18668, 5228, 18728, 548,
                                                                       563, 6008, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20018, 0, 3,
                                                                       18728, 5258, 18788, 563,
                                                                       578, 6053, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20108, 0, 3,
                                                                       18788, 5288, 18848, 578,
                                                                       593, 6098, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20198, 0, 3,
                                                                       18848, 5318, 18908, 593,
                                                                       608, 6143, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20288, 0, 3,
                                                                       18908, 5348, 18968, 608,
                                                                       623, 6188, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20378, 0, 3,
                                                                       18968, 5378, 19028, 623,
                                                                       638, 6233, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20468, 0, 3,
                                                                       19028, 5408, 19088, 638,
                                                                       653, 6278, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20558, 0, 3,
                                                                       19148, 5528, 19208, 683,
                                                                       698, 6413, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20648, 0, 3,
                                                                       19208, 5558, 19268, 698,
                                                                       713, 6458, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20738, 0, 3,
                                                                       19268, 5588, 19328, 713,
                                                                       728, 6503, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20828, 0, 3,
                                                                       19328, 5618, 19388, 728,
                                                                       743, 6548, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20918, 0, 3,
                                                                       19388, 5648, 19448, 743,
                                                                       758, 6593, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21008, 0, 3,
                                                                       19448, 5678, 19508, 758,
                                                                       773, 6638, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21098, 0, 3,
                                                                       19508, 5708, 19568, 773,
                                                                       788, 6683, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21188, 0, 3,
                                                                       19568, 5738, 19628, 788,
                                                                       803, 6728, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21278, 0, 3,
                                                                       19628, 5768, 19688, 803,
                                                                       818, 6773, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21368, 0, 3,
                                                                       19748, 5918, 19838, 848,
                                                                       869, 6944, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21494, 0, 3,
                                                                       19838, 5963, 19928, 869,
                                                                       890, 7007, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21620, 0, 3,
                                                                       19928, 6008, 20018, 890,
                                                                       911, 7070, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21746, 0, 3,
                                                                       20018, 6053, 20108, 911,
                                                                       932, 7133, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21872, 0, 3,
                                                                       20108, 6098, 20198, 932,
                                                                       953, 7196, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21998, 0, 3,
                                                                       20198, 6143, 20288, 953,
                                                                       974, 7259, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22124, 0, 3,
                                                                       20288, 6188, 20378, 974,
                                                                       995, 7322, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22250, 0, 3,
                                                                       20378, 6233, 20468, 995,
                                                                       1016, 7385, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22376, 0, 3,
                                                                       20558, 6413, 20648, 1058,
                                                                       1079, 7574, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22502, 0, 3,
                                                                       20648, 6458, 20738, 1079,
                                                                       1100, 7637, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22628, 0, 3,
                                                                       20738, 6503, 20828, 1100,
                                                                       1121, 7700, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22754, 0, 3,
                                                                       20828, 6548, 20918, 1121,
                                                                       1142, 7763, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22880, 0, 3,
                                                                       20918, 6593, 21008, 1142,
                                                                       1163, 7826, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23006, 0, 3,
                                                                       21008, 6638, 21098, 1163,
                                                                       1184, 7889, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23132, 0, 3,
                                                                       21098, 6683, 21188, 1184,
                                                                       1205, 7952, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23258, 0, 3,
                                                                       21188, 6728, 21278, 1205,
                                                                       1226, 8015, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23384, 0, 3,
                                                                       21368, 6944, 21494, 1268,
                                                                       1296, 8246, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23552, 0, 3,
                                                                       21494, 7007, 21620, 1296,
                                                                       1324, 8330, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23720, 0, 3,
                                                                       21620, 7070, 21746, 1324,
                                                                       1352, 8414, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23888, 0, 3,
                                                                       21746, 7133, 21872, 1352,
                                                                       1380, 8498, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24056, 0, 3,
                                                                       21872, 7196, 21998, 1380,
                                                                       1408, 8582, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24224, 0, 3,
                                                                       21998, 7259, 22124, 1408,
                                                                       1436, 8666, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24392, 0, 3,
                                                                       22124, 7322, 22250, 1436,
                                                                       1464, 8750, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24560, 0, 3,
                                                                       22376, 7574, 22502, 1520,
                                                                       1548, 9002, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24728, 0, 3,
                                                                       22502, 7637, 22628, 1548,
                                                                       1576, 9086, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24896, 0, 3,
                                                                       22628, 7700, 22754, 1576,
                                                                       1604, 9170, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25064, 0, 3,
                                                                       22754, 7763, 22880, 1604,
                                                                       1632, 9254, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25232, 0, 3,
                                                                       22880, 7826, 23006, 1632,
                                                                       1660, 9338, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25400, 0, 3,
                                                                       23006, 7889, 23132, 1660,
                                                                       1688, 9422, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25568, 0, 3,
                                                                       23132, 7952, 23258, 1688,
                                                                       1716, 9506, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25736, 0, 3,
                                                                       23384, 8246, 23552, 1772,
                                                                       1808, 9806, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25952, 0, 3,
                                                                       23552, 8330, 23720, 1808,
                                                                       1844, 9914, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26168, 0, 3,
                                                                       23720, 8414, 23888, 1844,
                                                                       1880, 10022, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26384, 0, 3,
                                                                       23888, 8498, 24056, 1880,
                                                                       1916, 10130, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26600, 0, 3,
                                                                       24056, 8582, 24224, 1916,
                                                                       1952, 10238, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26816, 0, 3,
                                                                       24224, 8666, 24392, 1952,
                                                                       1988, 10346, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27032, 0, 3,
                                                                       24560, 9002, 24728, 2060,
                                                                       2096, 10670, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27248, 0, 3,
                                                                       24728, 9086, 24896, 2096,
                                                                       2132, 10778, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27464, 0, 3,
                                                                       24896, 9170, 25064, 2132,
                                                                       2168, 10886, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27680, 0, 3,
                                                                       25064, 9254, 25232, 2168,
                                                                       2204, 10994, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27896, 0, 3,
                                                                       25232, 9338, 25400, 2204,
                                                                       2240, 11102, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28112, 0, 3,
                                                                       25400, 9422, 25568, 2240,
                                                                       2276, 11210, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28328, 0, 3,
                                                                       25736, 9806, 25952, 2348,
                                                                       2393, 11588, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28598, 0, 3,
                                                                       25952, 9914, 26168, 2393,
                                                                       2438, 11723, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28868, 0, 3,
                                                                       26168, 10022, 26384, 2438,
                                                                       2483, 11858, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29138, 0, 3,
                                                                       26384, 10130, 26600, 2483,
                                                                       2528, 11993, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29408, 0, 3,
                                                                       26600, 10238, 26816, 2528,
                                                                       2573, 12128, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29678, 0, 3,
                                                                       27032, 10670, 27248, 2663,
                                                                       2708, 12533, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29948, 0, 3,
                                                                       27248, 10778, 27464, 2708,
                                                                       2753, 12668, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 30218, 0, 3,
                                                                       27464, 10886, 27680, 2753,
                                                                       2798, 12803, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 30488, 0, 3,
                                                                       27680, 10994, 27896, 2798,
                                                                       2843, 12938, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 30758, 0, 3,
                                                                       27896, 11102, 28112, 2843,
                                                                       2888, 13073, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 31028, 0, 3,
                                                                       28328, 11588, 28598, 2978,
                                                                       3033, 13538, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 31358, 0, 3,
                                                                       28598, 11723, 28868, 3033,
                                                                       3088, 13703, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 31688, 0, 3,
                                                                       28868, 11858, 29138, 3088,
                                                                       3143, 13868, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 32018, 0, 3,
                                                                       29138, 11993, 29408, 3143,
                                                                       3198, 14033, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 32348, 0, 3,
                                                                       29678, 12533, 29948, 3308,
                                                                       3363, 14528, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 32678, 0, 3,
                                                                       29948, 12668, 30218, 3363,
                                                                       3418, 14693, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 33008, 0, 3,
                                                                       30218, 12803, 30488, 3418,
                                                                       3473, 14858, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 33338, 0, 3,
                                                                       30488, 12938, 30758, 3473,
                                                                       3528, 15023, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 33668, 0, 3,
                                                                       31028, 13538, 31358, 3638,
                                                                       3704, 15584, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 34064, 0, 3,
                                                                       31358, 13703, 31688, 3704,
                                                                       3770, 15782, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 34460, 0, 3,
                                                                       31688, 13868, 32018, 3770,
                                                                       3836, 15980, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 34856, 0, 3,
                                                                       32348, 14528, 32678, 3968,
                                                                       4034, 16574, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 35252, 0, 3,
                                                                       32678, 14693, 33008, 4034,
                                                                       4100, 16772, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 35648, 0, 3,
                                                                       33008, 14858, 33338, 4100,
                                                                       4166, 16970, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36044, 3, 4298,
                                                                       4301, 17168, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36054, 3, 4301,
                                                                       4304, 17174, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36064, 3, 4304,
                                                                       4307, 17180, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36074, 3, 4307,
                                                                       4310, 17186, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36084, 3, 4310,
                                                                       4313, 17192, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36094, 3, 4313,
                                                                       4316, 17198, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36104, 3, 4316,
                                                                       4319, 17204, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36114, 3, 4319,
                                                                       4322, 17210, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36124, 3, 4322,
                                                                       4325, 17216, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36134, 3, 4325,
                                                                       4328, 17222, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36144, 3, 4328,
                                                                       4331, 17228, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36154, 3, 4331,
                                                                       4334, 17234, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36164, 3, 4334,
                                                                       4337, 17240, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36174, 3, 4343,
                                                                       4346, 17246, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36184, 3, 4346,
                                                                       4349, 17252, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36194, 3, 4349,
                                                                       4352, 17258, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36204, 3, 4352,
                                                                       4355, 17264, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36214, 3, 4355,
                                                                       4358, 17270, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36224, 3, 4358,
                                                                       4361, 17276, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36234, 3, 4361,
                                                                       4364, 17282, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36244, 3, 4364,
                                                                       4367, 17288, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36254, 3, 4367,
                                                                       4370, 17294, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36264, 3, 4370,
                                                                       4373, 17300, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36274, 3, 4373,
                                                                       4376, 17306, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36284, 3, 4376,
                                                                       4379, 17312, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36294, 3, 4379,
                                                                       4382, 17318, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36304, 0, 3,
                                                                       36044, 17168, 36054, 4388,
                                                                       4397, 17324, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36334, 0, 3,
                                                                       36054, 17174, 36064, 4397,
                                                                       4406, 17342, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36364, 0, 3,
                                                                       36064, 17180, 36074, 4406,
                                                                       4415, 17360, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36394, 0, 3,
                                                                       36074, 17186, 36084, 4415,
                                                                       4424, 17378, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36424, 0, 3,
                                                                       36084, 17192, 36094, 4424,
                                                                       4433, 17396, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36454, 0, 3,
                                                                       36094, 17198, 36104, 4433,
                                                                       4442, 17414, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36484, 0, 3,
                                                                       36104, 17204, 36114, 4442,
                                                                       4451, 17432, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36514, 0, 3,
                                                                       36114, 17210, 36124, 4451,
                                                                       4460, 17450, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36544, 0, 3,
                                                                       36124, 17216, 36134, 4460,
                                                                       4469, 17468, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36574, 0, 3,
                                                                       36134, 17222, 36144, 4469,
                                                                       4478, 17486, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36604, 0, 3,
                                                                       36144, 17228, 36154, 4478,
                                                                       4487, 17504, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36634, 0, 3,
                                                                       36154, 17234, 36164, 4487,
                                                                       4496, 17522, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36664, 0, 3,
                                                                       36174, 17246, 36184, 4514,
                                                                       4523, 17540, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36694, 0, 3,
                                                                       36184, 17252, 36194, 4523,
                                                                       4532, 17558, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36724, 0, 3,
                                                                       36194, 17258, 36204, 4532,
                                                                       4541, 17576, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36754, 0, 3,
                                                                       36204, 17264, 36214, 4541,
                                                                       4550, 17594, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36784, 0, 3,
                                                                       36214, 17270, 36224, 4550,
                                                                       4559, 17612, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36814, 0, 3,
                                                                       36224, 17276, 36234, 4559,
                                                                       4568, 17630, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36844, 0, 3,
                                                                       36234, 17282, 36244, 4568,
                                                                       4577, 17648, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36874, 0, 3,
                                                                       36244, 17288, 36254, 4577,
                                                                       4586, 17666, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36904, 0, 3,
                                                                       36254, 17294, 36264, 4586,
                                                                       4595, 17684, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36934, 0, 3,
                                                                       36264, 17300, 36274, 4595,
                                                                       4604, 17702, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36964, 0, 3,
                                                                       36274, 17306, 36284, 4604,
                                                                       4613, 17720, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 36994, 0, 3,
                                                                       36284, 17312, 36294, 4613,
                                                                       4622, 17738, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37024, 0, 3,
                                                                       36304, 17324, 36334, 4640,
                                                                       4658, 17756, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37084, 0, 3,
                                                                       36334, 17342, 36364, 4658,
                                                                       4676, 17792, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37144, 0, 3,
                                                                       36364, 17360, 36394, 4676,
                                                                       4694, 17828, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37204, 0, 3,
                                                                       36394, 17378, 36424, 4694,
                                                                       4712, 17864, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37264, 0, 3,
                                                                       36424, 17396, 36454, 4712,
                                                                       4730, 17900, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37324, 0, 3,
                                                                       36454, 17414, 36484, 4730,
                                                                       4748, 17936, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37384, 0, 3,
                                                                       36484, 17432, 36514, 4748,
                                                                       4766, 17972, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37444, 0, 3,
                                                                       36514, 17450, 36544, 4766,
                                                                       4784, 18008, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37504, 0, 3,
                                                                       36544, 17468, 36574, 4784,
                                                                       4802, 18044, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37564, 0, 3,
                                                                       36574, 17486, 36604, 4802,
                                                                       4820, 18080, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37624, 0, 3,
                                                                       36604, 17504, 36634, 4820,
                                                                       4838, 18116, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37684, 0, 3,
                                                                       36664, 17540, 36694, 4874,
                                                                       4892, 18152, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37744, 0, 3,
                                                                       36694, 17558, 36724, 4892,
                                                                       4910, 18188, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37804, 0, 3,
                                                                       36724, 17576, 36754, 4910,
                                                                       4928, 18224, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37864, 0, 3,
                                                                       36754, 17594, 36784, 4928,
                                                                       4946, 18260, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37924, 0, 3,
                                                                       36784, 17612, 36814, 4946,
                                                                       4964, 18296, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37984, 0, 3,
                                                                       36814, 17630, 36844, 4964,
                                                                       4982, 18332, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38044, 0, 3,
                                                                       36844, 17648, 36874, 4982,
                                                                       5000, 18368, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38104, 0, 3,
                                                                       36874, 17666, 36904, 5000,
                                                                       5018, 18404, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38164, 0, 3,
                                                                       36904, 17684, 36934, 5018,
                                                                       5036, 18440, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38224, 0, 3,
                                                                       36934, 17702, 36964, 5036,
                                                                       5054, 18476, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38284, 0, 3,
                                                                       36964, 17720, 36994, 5054,
                                                                       5072, 18512, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38344, 0, 3,
                                                                       37024, 17756, 37084, 5108,
                                                                       5138, 18548, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38444, 0, 3,
                                                                       37084, 17792, 37144, 5138,
                                                                       5168, 18608, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38544, 0, 3,
                                                                       37144, 17828, 37204, 5168,
                                                                       5198, 18668, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38644, 0, 3,
                                                                       37204, 17864, 37264, 5198,
                                                                       5228, 18728, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38744, 0, 3,
                                                                       37264, 17900, 37324, 5228,
                                                                       5258, 18788, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38844, 0, 3,
                                                                       37324, 17936, 37384, 5258,
                                                                       5288, 18848, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38944, 0, 3,
                                                                       37384, 17972, 37444, 5288,
                                                                       5318, 18908, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39044, 0, 3,
                                                                       37444, 18008, 37504, 5318,
                                                                       5348, 18968, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39144, 0, 3,
                                                                       37504, 18044, 37564, 5348,
                                                                       5378, 19028, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39244, 0, 3,
                                                                       37564, 18080, 37624, 5378,
                                                                       5408, 19088, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39344, 0, 3,
                                                                       37684, 18152, 37744, 5468,
                                                                       5498, 19148, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39444, 0, 3,
                                                                       37744, 18188, 37804, 5498,
                                                                       5528, 19208, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39544, 0, 3,
                                                                       37804, 18224, 37864, 5528,
                                                                       5558, 19268, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39644, 0, 3,
                                                                       37864, 18260, 37924, 5558,
                                                                       5588, 19328, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39744, 0, 3,
                                                                       37924, 18296, 37984, 5588,
                                                                       5618, 19388, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39844, 0, 3,
                                                                       37984, 18332, 38044, 5618,
                                                                       5648, 19448, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39944, 0, 3,
                                                                       38044, 18368, 38104, 5648,
                                                                       5678, 19508, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 40044, 0, 3,
                                                                       38104, 18404, 38164, 5678,
                                                                       5708, 19568, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 40144, 0, 3,
                                                                       38164, 18440, 38224, 5708,
                                                                       5738, 19628, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 40244, 0, 3,
                                                                       38224, 18476, 38284, 5738,
                                                                       5768, 19688, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40344, 0, 3,
                                                                       38344, 18548, 38444, 5828,
                                                                       5873, 19748, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40494, 0, 3,
                                                                       38444, 18608, 38544, 5873,
                                                                       5918, 19838, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40644, 0, 3,
                                                                       38544, 18668, 38644, 5918,
                                                                       5963, 19928, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40794, 0, 3,
                                                                       38644, 18728, 38744, 5963,
                                                                       6008, 20018, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40944, 0, 3,
                                                                       38744, 18788, 38844, 6008,
                                                                       6053, 20108, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 41094, 0, 3,
                                                                       38844, 18848, 38944, 6053,
                                                                       6098, 20198, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 41244, 0, 3,
                                                                       38944, 18908, 39044, 6098,
                                                                       6143, 20288, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 41394, 0, 3,
                                                                       39044, 18968, 39144, 6143,
                                                                       6188, 20378, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 41544, 0, 3,
                                                                       39144, 19028, 39244, 6188,
                                                                       6233, 20468, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 41694, 0, 3,
                                                                       39344, 19148, 39444, 6323,
                                                                       6368, 20558, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 41844, 0, 3,
                                                                       39444, 19208, 39544, 6368,
                                                                       6413, 20648, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 41994, 0, 3,
                                                                       39544, 19268, 39644, 6413,
                                                                       6458, 20738, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 42144, 0, 3,
                                                                       39644, 19328, 39744, 6458,
                                                                       6503, 20828, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 42294, 0, 3,
                                                                       39744, 19388, 39844, 6503,
                                                                       6548, 20918, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 42444, 0, 3,
                                                                       39844, 19448, 39944, 6548,
                                                                       6593, 21008, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 42594, 0, 3,
                                                                       39944, 19508, 40044, 6593,
                                                                       6638, 21098, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 42744, 0, 3,
                                                                       40044, 19568, 40144, 6638,
                                                                       6683, 21188, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 42894, 0, 3,
                                                                       40144, 19628, 40244, 6683,
                                                                       6728, 21278, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 43044, 0, 3,
                                                                       40344, 19748, 40494, 6818,
                                                                       6881, 21368, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 43254, 0, 3,
                                                                       40494, 19838, 40644, 6881,
                                                                       6944, 21494, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 43464, 0, 3,
                                                                       40644, 19928, 40794, 6944,
                                                                       7007, 21620, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 43674, 0, 3,
                                                                       40794, 20018, 40944, 7007,
                                                                       7070, 21746, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 43884, 0, 3,
                                                                       40944, 20108, 41094, 7070,
                                                                       7133, 21872, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 44094, 0, 3,
                                                                       41094, 20198, 41244, 7133,
                                                                       7196, 21998, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 44304, 0, 3,
                                                                       41244, 20288, 41394, 7196,
                                                                       7259, 22124, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 44514, 0, 3,
                                                                       41394, 20378, 41544, 7259,
                                                                       7322, 22250, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 44724, 0, 3,
                                                                       41694, 20558, 41844, 7448,
                                                                       7511, 22376, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 44934, 0, 3,
                                                                       41844, 20648, 41994, 7511,
                                                                       7574, 22502, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 45144, 0, 3,
                                                                       41994, 20738, 42144, 7574,
                                                                       7637, 22628, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 45354, 0, 3,
                                                                       42144, 20828, 42294, 7637,
                                                                       7700, 22754, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 45564, 0, 3,
                                                                       42294, 20918, 42444, 7700,
                                                                       7763, 22880, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 45774, 0, 3,
                                                                       42444, 21008, 42594, 7763,
                                                                       7826, 23006, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 45984, 0, 3,
                                                                       42594, 21098, 42744, 7826,
                                                                       7889, 23132, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 46194, 0, 3,
                                                                       42744, 21188, 42894, 7889,
                                                                       7952, 23258, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 46404, 0, 3,
                                                                       43044, 21368, 43254, 8078,
                                                                       8162, 23384, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 46684, 0, 3,
                                                                       43254, 21494, 43464, 8162,
                                                                       8246, 23552, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 46964, 0, 3,
                                                                       43464, 21620, 43674, 8246,
                                                                       8330, 23720, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 47244, 0, 3,
                                                                       43674, 21746, 43884, 8330,
                                                                       8414, 23888, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 47524, 0, 3,
                                                                       43884, 21872, 44094, 8414,
                                                                       8498, 24056, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 47804, 0, 3,
                                                                       44094, 21998, 44304, 8498,
                                                                       8582, 24224, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 48084, 0, 3,
                                                                       44304, 22124, 44514, 8582,
                                                                       8666, 24392, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 48364, 0, 3,
                                                                       44724, 22376, 44934, 8834,
                                                                       8918, 24560, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 48644, 0, 3,
                                                                       44934, 22502, 45144, 8918,
                                                                       9002, 24728, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 48924, 0, 3,
                                                                       45144, 22628, 45354, 9002,
                                                                       9086, 24896, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 49204, 0, 3,
                                                                       45354, 22754, 45564, 9086,
                                                                       9170, 25064, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 49484, 0, 3,
                                                                       45564, 22880, 45774, 9170,
                                                                       9254, 25232, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 49764, 0, 3,
                                                                       45774, 23006, 45984, 9254,
                                                                       9338, 25400, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 50044, 0, 3,
                                                                       45984, 23132, 46194, 9338,
                                                                       9422, 25568, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 50324, 0, 3,
                                                                       46404, 23384, 46684, 9590,
                                                                       9698, 25736, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 50684, 0, 3,
                                                                       46684, 23552, 46964, 9698,
                                                                       9806, 25952, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 51044, 0, 3,
                                                                       46964, 23720, 47244, 9806,
                                                                       9914, 26168, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 51404, 0, 3,
                                                                       47244, 23888, 47524, 9914,
                                                                       10022, 26384, ncols,
                                                                       gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 51764, 0, 3,
                                                                       47524, 24056, 47804,
                                                                       10022, 10130, 26600,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 52124, 0, 3,
                                                                       47804, 24224, 48084,
                                                                       10130, 10238, 26816,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 52484, 0, 3,
                                                                       48364, 24560, 48644,
                                                                       10454, 10562, 27032,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 52844, 0, 3,
                                                                       48644, 24728, 48924,
                                                                       10562, 10670, 27248,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 53204, 0, 3,
                                                                       48924, 24896, 49204,
                                                                       10670, 10778, 27464,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 53564, 0, 3,
                                                                       49204, 25064, 49484,
                                                                       10778, 10886, 27680,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 53924, 0, 3,
                                                                       49484, 25232, 49764,
                                                                       10886, 10994, 27896,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 54284, 0, 3,
                                                                       49764, 25400, 50044,
                                                                       10994, 11102, 28112,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 54644, 0, 3,
                                                                       50324, 25736, 50684,
                                                                       11318, 11453, 28328,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 55094, 0, 3,
                                                                       50684, 25952, 51044,
                                                                       11453, 11588, 28598,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 55544, 0, 3,
                                                                       51044, 26168, 51404,
                                                                       11588, 11723, 28868,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 55994, 0, 3,
                                                                       51404, 26384, 51764,
                                                                       11723, 11858, 29138,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 56444, 0, 3,
                                                                       51764, 26600, 52124,
                                                                       11858, 11993, 29408,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 56894, 0, 3,
                                                                       52484, 27032, 52844,
                                                                       12263, 12398, 29678,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 57344, 0, 3,
                                                                       52844, 27248, 53204,
                                                                       12398, 12533, 29948,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 57794, 0, 3,
                                                                       53204, 27464, 53564,
                                                                       12533, 12668, 30218,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 58244, 0, 3,
                                                                       53564, 27680, 53924,
                                                                       12668, 12803, 30488,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 58694, 0, 3,
                                                                       53924, 27896, 54284,
                                                                       12803, 12938, 30758,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 59144, 0, 3,
                                                                       54644, 28328, 55094,
                                                                       13208, 13373, 31028,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 59694, 0, 3,
                                                                       55094, 28598, 55544,
                                                                       13373, 13538, 31358,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 60244, 0, 3,
                                                                       55544, 28868, 55994,
                                                                       13538, 13703, 31688,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 60794, 0, 3,
                                                                       55994, 29138, 56444,
                                                                       13703, 13868, 32018,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 61344, 0, 3,
                                                                       56894, 29678, 57344,
                                                                       14198, 14363, 32348,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 61894, 0, 3,
                                                                       57344, 29948, 57794,
                                                                       14363, 14528, 32678,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 62444, 0, 3,
                                                                       57794, 30218, 58244,
                                                                       14528, 14693, 33008,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 62994, 0, 3,
                                                                       58244, 30488, 58694,
                                                                       14693, 14858, 33338,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 63544, 0, 3,
                                                                       59144, 31028, 59694,
                                                                       15188, 15386, 33668,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 64204, 0, 3,
                                                                       59694, 31358, 60244,
                                                                       15386, 15584, 34064,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 64864, 0, 3,
                                                                       60244, 31688, 60794,
                                                                       15584, 15782, 34460,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 65524, 0, 3,
                                                                       61344, 32348, 61894,
                                                                       16178, 16376, 34856,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 66184, 0, 3,
                                                                       61894, 32678, 62444,
                                                                       16376, 16574, 35252,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 66844, 0, 3,
                                                                       62444, 33008, 62994,
                                                                       16574, 16772, 35648,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67504, 3, 17168,
                                                                       17174, 36064, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67519, 3, 17174,
                                                                       17180, 36074, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67534, 3, 17180,
                                                                       17186, 36084, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67549, 3, 17186,
                                                                       17192, 36094, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67564, 3, 17192,
                                                                       17198, 36104, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67579, 3, 17198,
                                                                       17204, 36114, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67594, 3, 17204,
                                                                       17210, 36124, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67609, 3, 17210,
                                                                       17216, 36134, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67624, 3, 17216,
                                                                       17222, 36144, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67639, 3, 17222,
                                                                       17228, 36154, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67654, 3, 17228,
                                                                       17234, 36164, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67669, 3, 17246,
                                                                       17252, 36194, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67684, 3, 17252,
                                                                       17258, 36204, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67699, 3, 17258,
                                                                       17264, 36214, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67714, 3, 17264,
                                                                       17270, 36224, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67729, 3, 17270,
                                                                       17276, 36234, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67744, 3, 17276,
                                                                       17282, 36244, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67759, 3, 17282,
                                                                       17288, 36254, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67774, 3, 17288,
                                                                       17294, 36264, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67789, 3, 17294,
                                                                       17300, 36274, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67804, 3, 17300,
                                                                       17306, 36284, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67819, 3, 17306,
                                                                       17312, 36294, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 67834, 0, 3,
                                                                       67504, 36064, 67519,
                                                                       17324, 17342, 36364,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 67879, 0, 3,
                                                                       67519, 36074, 67534,
                                                                       17342, 17360, 36394,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 67924, 0, 3,
                                                                       67534, 36084, 67549,
                                                                       17360, 17378, 36424,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 67969, 0, 3,
                                                                       67549, 36094, 67564,
                                                                       17378, 17396, 36454,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68014, 0, 3,
                                                                       67564, 36104, 67579,
                                                                       17396, 17414, 36484,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68059, 0, 3,
                                                                       67579, 36114, 67594,
                                                                       17414, 17432, 36514,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68104, 0, 3,
                                                                       67594, 36124, 67609,
                                                                       17432, 17450, 36544,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68149, 0, 3,
                                                                       67609, 36134, 67624,
                                                                       17450, 17468, 36574,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68194, 0, 3,
                                                                       67624, 36144, 67639,
                                                                       17468, 17486, 36604,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68239, 0, 3,
                                                                       67639, 36154, 67654,
                                                                       17486, 17504, 36634,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68284, 0, 3,
                                                                       67669, 36194, 67684,
                                                                       17540, 17558, 36724,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68329, 0, 3,
                                                                       67684, 36204, 67699,
                                                                       17558, 17576, 36754,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68374, 0, 3,
                                                                       67699, 36214, 67714,
                                                                       17576, 17594, 36784,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68419, 0, 3,
                                                                       67714, 36224, 67729,
                                                                       17594, 17612, 36814,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68464, 0, 3,
                                                                       67729, 36234, 67744,
                                                                       17612, 17630, 36844,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68509, 0, 3,
                                                                       67744, 36244, 67759,
                                                                       17630, 17648, 36874,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68554, 0, 3,
                                                                       67759, 36254, 67774,
                                                                       17648, 17666, 36904,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68599, 0, 3,
                                                                       67774, 36264, 67789,
                                                                       17666, 17684, 36934,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68644, 0, 3,
                                                                       67789, 36274, 67804,
                                                                       17684, 17702, 36964,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 68689, 0, 3,
                                                                       67804, 36284, 67819,
                                                                       17702, 17720, 36994,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 68734, 0, 3,
                                                                       67834, 36364, 67879,
                                                                       17756, 17792, 37144,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 68824, 0, 3,
                                                                       67879, 36394, 67924,
                                                                       17792, 17828, 37204,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 68914, 0, 3,
                                                                       67924, 36424, 67969,
                                                                       17828, 17864, 37264,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69004, 0, 3,
                                                                       67969, 36454, 68014,
                                                                       17864, 17900, 37324,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69094, 0, 3,
                                                                       68014, 36484, 68059,
                                                                       17900, 17936, 37384,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69184, 0, 3,
                                                                       68059, 36514, 68104,
                                                                       17936, 17972, 37444,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69274, 0, 3,
                                                                       68104, 36544, 68149,
                                                                       17972, 18008, 37504,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69364, 0, 3,
                                                                       68149, 36574, 68194,
                                                                       18008, 18044, 37564,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69454, 0, 3,
                                                                       68194, 36604, 68239,
                                                                       18044, 18080, 37624,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69544, 0, 3,
                                                                       68284, 36724, 68329,
                                                                       18152, 18188, 37804,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69634, 0, 3,
                                                                       68329, 36754, 68374,
                                                                       18188, 18224, 37864,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69724, 0, 3,
                                                                       68374, 36784, 68419,
                                                                       18224, 18260, 37924,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69814, 0, 3,
                                                                       68419, 36814, 68464,
                                                                       18260, 18296, 37984,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69904, 0, 3,
                                                                       68464, 36844, 68509,
                                                                       18296, 18332, 38044,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 69994, 0, 3,
                                                                       68509, 36874, 68554,
                                                                       18332, 18368, 38104,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 70084, 0, 3,
                                                                       68554, 36904, 68599,
                                                                       18368, 18404, 38164,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 70174, 0, 3,
                                                                       68599, 36934, 68644,
                                                                       18404, 18440, 38224,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 70264, 0, 3,
                                                                       68644, 36964, 68689,
                                                                       18440, 18476, 38284,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 70354, 0, 3,
                                                                       68734, 37144, 68824,
                                                                       18548, 18608, 38544,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 70504, 0, 3,
                                                                       68824, 37204, 68914,
                                                                       18608, 18668, 38644,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 70654, 0, 3,
                                                                       68914, 37264, 69004,
                                                                       18668, 18728, 38744,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 70804, 0, 3,
                                                                       69004, 37324, 69094,
                                                                       18728, 18788, 38844,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 70954, 0, 3,
                                                                       69094, 37384, 69184,
                                                                       18788, 18848, 38944,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 71104, 0, 3,
                                                                       69184, 37444, 69274,
                                                                       18848, 18908, 39044,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 71254, 0, 3,
                                                                       69274, 37504, 69364,
                                                                       18908, 18968, 39144,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 71404, 0, 3,
                                                                       69364, 37564, 69454,
                                                                       18968, 19028, 39244,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 71554, 0, 3,
                                                                       69544, 37804, 69634,
                                                                       19148, 19208, 39544,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 71704, 0, 3,
                                                                       69634, 37864, 69724,
                                                                       19208, 19268, 39644,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 71854, 0, 3,
                                                                       69724, 37924, 69814,
                                                                       19268, 19328, 39744,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 72004, 0, 3,
                                                                       69814, 37984, 69904,
                                                                       19328, 19388, 39844,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 72154, 0, 3,
                                                                       69904, 38044, 69994,
                                                                       19388, 19448, 39944,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 72304, 0, 3,
                                                                       69994, 38104, 70084,
                                                                       19448, 19508, 40044,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 72454, 0, 3,
                                                                       70084, 38164, 70174,
                                                                       19508, 19568, 40144,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 72604, 0, 3,
                                                                       70174, 38224, 70264,
                                                                       19568, 19628, 40244,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 72754, 0, 3,
                                                                       70354, 38544, 70504,
                                                                       19748, 19838, 40644,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 72979, 0, 3,
                                                                       70504, 38644, 70654,
                                                                       19838, 19928, 40794,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 73204, 0, 3,
                                                                       70654, 38744, 70804,
                                                                       19928, 20018, 40944,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 73429, 0, 3,
                                                                       70804, 38844, 70954,
                                                                       20018, 20108, 41094,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 73654, 0, 3,
                                                                       70954, 38944, 71104,
                                                                       20108, 20198, 41244,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 73879, 0, 3,
                                                                       71104, 39044, 71254,
                                                                       20198, 20288, 41394,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 74104, 0, 3,
                                                                       71254, 39144, 71404,
                                                                       20288, 20378, 41544,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 74329, 0, 3,
                                                                       71554, 39544, 71704,
                                                                       20558, 20648, 41994,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 74554, 0, 3,
                                                                       71704, 39644, 71854,
                                                                       20648, 20738, 42144,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 74779, 0, 3,
                                                                       71854, 39744, 72004,
                                                                       20738, 20828, 42294,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 75004, 0, 3,
                                                                       72004, 39844, 72154,
                                                                       20828, 20918, 42444,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 75229, 0, 3,
                                                                       72154, 39944, 72304,
                                                                       20918, 21008, 42594,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 75454, 0, 3,
                                                                       72304, 40044, 72454,
                                                                       21008, 21098, 42744,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 75679, 0, 3,
                                                                       72454, 40144, 72604,
                                                                       21098, 21188, 42894,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 75904, 0, 3,
                                                                       72754, 40644, 72979,
                                                                       21368, 21494, 43464,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 76219, 0, 3,
                                                                       72979, 40794, 73204,
                                                                       21494, 21620, 43674,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 76534, 0, 3,
                                                                       73204, 40944, 73429,
                                                                       21620, 21746, 43884,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 76849, 0, 3,
                                                                       73429, 41094, 73654,
                                                                       21746, 21872, 44094,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 77164, 0, 3,
                                                                       73654, 41244, 73879,
                                                                       21872, 21998, 44304,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 77479, 0, 3,
                                                                       73879, 41394, 74104,
                                                                       21998, 22124, 44514,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 77794, 0, 3,
                                                                       74329, 41994, 74554,
                                                                       22376, 22502, 45144,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 78109, 0, 3,
                                                                       74554, 42144, 74779,
                                                                       22502, 22628, 45354,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 78424, 0, 3,
                                                                       74779, 42294, 75004,
                                                                       22628, 22754, 45564,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 78739, 0, 3,
                                                                       75004, 42444, 75229,
                                                                       22754, 22880, 45774,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 79054, 0, 3,
                                                                       75229, 42594, 75454,
                                                                       22880, 23006, 45984,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 79369, 0, 3,
                                                                       75454, 42744, 75679,
                                                                       23006, 23132, 46194,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 79684, 0, 3,
                                                                       75904, 43464, 76219,
                                                                       23384, 23552, 46964,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 80104, 0, 3,
                                                                       76219, 43674, 76534,
                                                                       23552, 23720, 47244,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 80524, 0, 3,
                                                                       76534, 43884, 76849,
                                                                       23720, 23888, 47524,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 80944, 0, 3,
                                                                       76849, 44094, 77164,
                                                                       23888, 24056, 47804,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 81364, 0, 3,
                                                                       77164, 44304, 77479,
                                                                       24056, 24224, 48084,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 81784, 0, 3,
                                                                       77794, 45144, 78109,
                                                                       24560, 24728, 48924,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 82204, 0, 3,
                                                                       78109, 45354, 78424,
                                                                       24728, 24896, 49204,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 82624, 0, 3,
                                                                       78424, 45564, 78739,
                                                                       24896, 25064, 49484,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 83044, 0, 3,
                                                                       78739, 45774, 79054,
                                                                       25064, 25232, 49764,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 83464, 0, 3,
                                                                       79054, 45984, 79369,
                                                                       25232, 25400, 50044,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 83884, 0, 3,
                                                                       79684, 46964, 80104,
                                                                       25736, 25952, 51044,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 84424, 0, 3,
                                                                       80104, 47244, 80524,
                                                                       25952, 26168, 51404,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 84964, 0, 3,
                                                                       80524, 47524, 80944,
                                                                       26168, 26384, 51764,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 85504, 0, 3,
                                                                       80944, 47804, 81364,
                                                                       26384, 26600, 52124,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 86044, 0, 3,
                                                                       81784, 48924, 82204,
                                                                       27032, 27248, 53204,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 86584, 0, 3,
                                                                       82204, 49204, 82624,
                                                                       27248, 27464, 53564,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 87124, 0, 3,
                                                                       82624, 49484, 83044,
                                                                       27464, 27680, 53924,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 87664, 0, 3,
                                                                       83044, 49764, 83464,
                                                                       27680, 27896, 54284,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 88204, 0, 3,
                                                                       83884, 51044, 84424,
                                                                       28328, 28598, 55544,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 88879, 0, 3,
                                                                       84424, 51404, 84964,
                                                                       28598, 28868, 55994,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 89554, 0, 3,
                                                                       84964, 51764, 85504,
                                                                       28868, 29138, 56444,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 90229, 0, 3,
                                                                       86044, 53204, 86584,
                                                                       29678, 29948, 57794,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 90904, 0, 3,
                                                                       86584, 53564, 87124,
                                                                       29948, 30218, 58244,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 91579, 0, 3,
                                                                       87124, 53924, 87664,
                                                                       30218, 30488, 58694,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 92254, 0, 3,
                                                                       88204, 55544, 88879,
                                                                       31028, 31358, 60244,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 93079, 0, 3,
                                                                       88879, 55994, 89554,
                                                                       31358, 31688, 60794,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 93904, 0, 3,
                                                                       90229, 57794, 90904,
                                                                       32348, 32678, 62444,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 94729, 0, 3,
                                                                       90904, 58244, 91579,
                                                                       32678, 33008, 62994,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 95554, 0, 3,
                                                                       92254, 60244, 93079,
                                                                       33668, 34064, 64864,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 96544, 0, 3,
                                                                       93904, 62444, 94729,
                                                                       34856, 35252, 66844,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97534, 3, 36044,
                                                                       36054, 67504, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97555, 3, 36054,
                                                                       36064, 67519, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97576, 3, 36064,
                                                                       36074, 67534, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97597, 3, 36074,
                                                                       36084, 67549, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97618, 3, 36084,
                                                                       36094, 67564, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97639, 3, 36094,
                                                                       36104, 67579, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97660, 3, 36104,
                                                                       36114, 67594, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97681, 3, 36114,
                                                                       36124, 67609, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97702, 3, 36124,
                                                                       36134, 67624, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97723, 3, 36134,
                                                                       36144, 67639, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97744, 3, 36144,
                                                                       36154, 67654, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97765, 3, 36174,
                                                                       36184, 67669, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97786, 3, 36184,
                                                                       36194, 67684, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97807, 3, 36194,
                                                                       36204, 67699, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97828, 3, 36204,
                                                                       36214, 67714, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97849, 3, 36214,
                                                                       36224, 67729, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97870, 3, 36224,
                                                                       36234, 67744, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97891, 3, 36234,
                                                                       36244, 67759, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97912, 3, 36244,
                                                                       36254, 67774, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97933, 3, 36254,
                                                                       36264, 67789, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97954, 3, 36264,
                                                                       36274, 67804, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 97975, 3, 36274,
                                                                       36284, 67819, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 97996, 0, 3,
                                                                       97534, 67504, 97555,
                                                                       36304, 36334, 67834,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 98059, 0, 3,
                                                                       97555, 67519, 97576,
                                                                       36334, 36364, 67879,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 98122, 0, 3,
                                                                       97576, 67534, 97597,
                                                                       36364, 36394, 67924,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 98185, 0, 3,
                                                                       97597, 67549, 97618,
                                                                       36394, 36424, 67969,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 98248, 0, 3,
                                                                       97618, 67564, 97639,
                                                                       36424, 36454, 68014,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 98311, 0, 3,
                                                                       97639, 67579, 97660,
                                                                       36454, 36484, 68059,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 98374, 0, 3,
                                                                       97660, 67594, 97681,
                                                                       36484, 36514, 68104,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 98437, 0, 3,
                                                                       97681, 67609, 97702,
                                                                       36514, 36544, 68149,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 98500, 0, 3,
                                                                       97702, 67624, 97723,
                                                                       36544, 36574, 68194,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 98563, 0, 3,
                                                                       97723, 67639, 97744,
                                                                       36574, 36604, 68239,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 98626, 0, 3,
                                                                       97765, 67669, 97786,
                                                                       36664, 36694, 68284,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 98689, 0, 3,
                                                                       97786, 67684, 97807,
                                                                       36694, 36724, 68329,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 98752, 0, 3,
                                                                       97807, 67699, 97828,
                                                                       36724, 36754, 68374,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 98815, 0, 3,
                                                                       97828, 67714, 97849,
                                                                       36754, 36784, 68419,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 98878, 0, 3,
                                                                       97849, 67729, 97870,
                                                                       36784, 36814, 68464,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 98941, 0, 3,
                                                                       97870, 67744, 97891,
                                                                       36814, 36844, 68509,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 99004, 0, 3,
                                                                       97891, 67759, 97912,
                                                                       36844, 36874, 68554,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 99067, 0, 3,
                                                                       97912, 67774, 97933,
                                                                       36874, 36904, 68599,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 99130, 0, 3,
                                                                       97933, 67789, 97954,
                                                                       36904, 36934, 68644,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 99193, 0, 3,
                                                                       97954, 67804, 97975,
                                                                       36934, 36964, 68689,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 99256, 0, 3,
                                                                       97996, 67834, 98059,
                                                                       37024, 37084, 68734,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 99382, 0, 3,
                                                                       98059, 67879, 98122,
                                                                       37084, 37144, 68824,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 99508, 0, 3,
                                                                       98122, 67924, 98185,
                                                                       37144, 37204, 68914,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 99634, 0, 3,
                                                                       98185, 67969, 98248,
                                                                       37204, 37264, 69004,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 99760, 0, 3,
                                                                       98248, 68014, 98311,
                                                                       37264, 37324, 69094,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 99886, 0, 3,
                                                                       98311, 68059, 98374,
                                                                       37324, 37384, 69184,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 100012, 0, 3,
                                                                       98374, 68104, 98437,
                                                                       37384, 37444, 69274,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 100138, 0, 3,
                                                                       98437, 68149, 98500,
                                                                       37444, 37504, 69364,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 100264, 0, 3,
                                                                       98500, 68194, 98563,
                                                                       37504, 37564, 69454,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 100390, 0, 3,
                                                                       98626, 68284, 98689,
                                                                       37684, 37744, 69544,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 100516, 0, 3,
                                                                       98689, 68329, 98752,
                                                                       37744, 37804, 69634,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 100642, 0, 3,
                                                                       98752, 68374, 98815,
                                                                       37804, 37864, 69724,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 100768, 0, 3,
                                                                       98815, 68419, 98878,
                                                                       37864, 37924, 69814,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 100894, 0, 3,
                                                                       98878, 68464, 98941,
                                                                       37924, 37984, 69904,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 101020, 0, 3,
                                                                       98941, 68509, 99004,
                                                                       37984, 38044, 69994,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 101146, 0, 3,
                                                                       99004, 68554, 99067,
                                                                       38044, 38104, 70084,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 101272, 0, 3,
                                                                       99067, 68599, 99130,
                                                                       38104, 38164, 70174,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 101398, 0, 3,
                                                                       99130, 68644, 99193,
                                                                       38164, 38224, 70264,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 101524, 0, 3,
                                                                       99256, 68734, 99382,
                                                                       38344, 38444, 70354,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 101734, 0, 3,
                                                                       99382, 68824, 99508,
                                                                       38444, 38544, 70504,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 101944, 0, 3,
                                                                       99508, 68914, 99634,
                                                                       38544, 38644, 70654,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 102154, 0, 3,
                                                                       99634, 69004, 99760,
                                                                       38644, 38744, 70804,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 102364, 0, 3,
                                                                       99760, 69094, 99886,
                                                                       38744, 38844, 70954,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 102574, 0, 3,
                                                                       99886, 69184, 100012,
                                                                       38844, 38944, 71104,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 102784, 0, 3,
                                                                       100012, 69274, 100138,
                                                                       38944, 39044, 71254,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 102994, 0, 3,
                                                                       100138, 69364, 100264,
                                                                       39044, 39144, 71404,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 103204, 0, 3,
                                                                       100390, 69544, 100516,
                                                                       39344, 39444, 71554,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 103414, 0, 3,
                                                                       100516, 69634, 100642,
                                                                       39444, 39544, 71704,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 103624, 0, 3,
                                                                       100642, 69724, 100768,
                                                                       39544, 39644, 71854,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 103834, 0, 3,
                                                                       100768, 69814, 100894,
                                                                       39644, 39744, 72004,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 104044, 0, 3,
                                                                       100894, 69904, 101020,
                                                                       39744, 39844, 72154,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 104254, 0, 3,
                                                                       101020, 69994, 101146,
                                                                       39844, 39944, 72304,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 104464, 0, 3,
                                                                       101146, 70084, 101272,
                                                                       39944, 40044, 72454,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 104674, 0, 3,
                                                                       101272, 70174, 101398,
                                                                       40044, 40144, 72604,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 104884, 0, 3,
                                                                       101524, 70354, 101734,
                                                                       40344, 40494, 72754,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 105199, 0, 3,
                                                                       101734, 70504, 101944,
                                                                       40494, 40644, 72979,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 105514, 0, 3,
                                                                       101944, 70654, 102154,
                                                                       40644, 40794, 73204,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 105829, 0, 3,
                                                                       102154, 70804, 102364,
                                                                       40794, 40944, 73429,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 106144, 0, 3,
                                                                       102364, 70954, 102574,
                                                                       40944, 41094, 73654,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 106459, 0, 3,
                                                                       102574, 71104, 102784,
                                                                       41094, 41244, 73879,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 106774, 0, 3,
                                                                       102784, 71254, 102994,
                                                                       41244, 41394, 74104,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 107089, 0, 3,
                                                                       103204, 71554, 103414,
                                                                       41694, 41844, 74329,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 107404, 0, 3,
                                                                       103414, 71704, 103624,
                                                                       41844, 41994, 74554,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 107719, 0, 3,
                                                                       103624, 71854, 103834,
                                                                       41994, 42144, 74779,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 108034, 0, 3,
                                                                       103834, 72004, 104044,
                                                                       42144, 42294, 75004,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 108349, 0, 3,
                                                                       104044, 72154, 104254,
                                                                       42294, 42444, 75229,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 108664, 0, 3,
                                                                       104254, 72304, 104464,
                                                                       42444, 42594, 75454,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 108979, 0, 3,
                                                                       104464, 72454, 104674,
                                                                       42594, 42744, 75679,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 109294, 0, 3,
                                                                       104884, 72754, 105199,
                                                                       43044, 43254, 75904,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 109735, 0, 3,
                                                                       105199, 72979, 105514,
                                                                       43254, 43464, 76219,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 110176, 0, 3,
                                                                       105514, 73204, 105829,
                                                                       43464, 43674, 76534,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 110617, 0, 3,
                                                                       105829, 73429, 106144,
                                                                       43674, 43884, 76849,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 111058, 0, 3,
                                                                       106144, 73654, 106459,
                                                                       43884, 44094, 77164,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 111499, 0, 3,
                                                                       106459, 73879, 106774,
                                                                       44094, 44304, 77479,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 111940, 0, 3,
                                                                       107089, 74329, 107404,
                                                                       44724, 44934, 77794,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 112381, 0, 3,
                                                                       107404, 74554, 107719,
                                                                       44934, 45144, 78109,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 112822, 0, 3,
                                                                       107719, 74779, 108034,
                                                                       45144, 45354, 78424,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 113263, 0, 3,
                                                                       108034, 75004, 108349,
                                                                       45354, 45564, 78739,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 113704, 0, 3,
                                                                       108349, 75229, 108664,
                                                                       45564, 45774, 79054,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 114145, 0, 3,
                                                                       108664, 75454, 108979,
                                                                       45774, 45984, 79369,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 114586, 0, 3,
                                                                       109294, 75904, 109735,
                                                                       46404, 46684, 79684,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 115174, 0, 3,
                                                                       109735, 76219, 110176,
                                                                       46684, 46964, 80104,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 115762, 0, 3,
                                                                       110176, 76534, 110617,
                                                                       46964, 47244, 80524,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 116350, 0, 3,
                                                                       110617, 76849, 111058,
                                                                       47244, 47524, 80944,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 116938, 0, 3,
                                                                       111058, 77164, 111499,
                                                                       47524, 47804, 81364,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 117526, 0, 3,
                                                                       111940, 77794, 112381,
                                                                       48364, 48644, 81784,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 118114, 0, 3,
                                                                       112381, 78109, 112822,
                                                                       48644, 48924, 82204,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 118702, 0, 3,
                                                                       112822, 78424, 113263,
                                                                       48924, 49204, 82624,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 119290, 0, 3,
                                                                       113263, 78739, 113704,
                                                                       49204, 49484, 83044,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 119878, 0, 3,
                                                                       113704, 79054, 114145,
                                                                       49484, 49764, 83464,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 120466, 0, 3,
                                                                       114586, 79684, 115174,
                                                                       50324, 50684, 83884,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 121222, 0, 3,
                                                                       115174, 80104, 115762,
                                                                       50684, 51044, 84424,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 121978, 0, 3,
                                                                       115762, 80524, 116350,
                                                                       51044, 51404, 84964,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 122734, 0, 3,
                                                                       116350, 80944, 116938,
                                                                       51404, 51764, 85504,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 123490, 0, 3,
                                                                       117526, 81784, 118114,
                                                                       52484, 52844, 86044,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 124246, 0, 3,
                                                                       118114, 82204, 118702,
                                                                       52844, 53204, 86584,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 125002, 0, 3,
                                                                       118702, 82624, 119290,
                                                                       53204, 53564, 87124,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 125758, 0, 3,
                                                                       119290, 83044, 119878,
                                                                       53564, 53924, 87664,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 126514, 0, 3,
                                                                       120466, 83884, 121222,
                                                                       54644, 55094, 88204,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 127459, 0, 3,
                                                                       121222, 84424, 121978,
                                                                       55094, 55544, 88879,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 128404, 0, 3,
                                                                       121978, 84964, 122734,
                                                                       55544, 55994, 89554,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 129349, 0, 3,
                                                                       123490, 86044, 124246,
                                                                       56894, 57344, 90229,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 130294, 0, 3,
                                                                       124246, 86584, 125002,
                                                                       57344, 57794, 90904,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 131239, 0, 3,
                                                                       125002, 87124, 125758,
                                                                       57794, 58244, 91579,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 132184, 0, 3,
                                                                       126514, 88204, 127459,
                                                                       59144, 59694, 92254,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 133339, 0, 3,
                                                                       127459, 88879, 128404,
                                                                       59694, 60244, 93079,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 134494, 0, 3,
                                                                       129349, 90229, 130294,
                                                                       61344, 61894, 93904,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 135649, 0, 3,
                                                                       130294, 90904, 131239,
                                                                       61894, 62444, 94729,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 136804, 0, 3,
                                                                       132184, 92254, 133339,
                                                                       63544, 64204, 95554,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 138190, 0, 3,
                                                                       134494, 93904, 135649,
                                                                       65524, 66184, 96544,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 139576, 109294, 441, ncols);

                    simdfunc::contract_primitives(buffer, 140248, 111940, 441, ncols);

                    simdfunc::contract_primitives(buffer, 140920, 114586, 588, ncols);

                    simdfunc::contract_primitives(buffer, 141816, 117526, 588, ncols);

                    simdfunc::contract_primitives(buffer, 142712, 120466, 756, ncols);

                    simdfunc::contract_primitives(buffer, 143864, 123490, 756, ncols);

                    simdfunc::contract_primitives(buffer, 145016, 126514, 945, ncols);

                    simdfunc::contract_primitives(buffer, 146456, 129349, 945, ncols);

                    simdfunc::contract_primitives(buffer, 147896, 132184, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 149656, 134494, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 151416, 136804, 1386, ncols);

                    simdfunc::contract_primitives(buffer, 153528, 138190, 1386, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 140017, 139576, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 140689, 140248, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 141508, 140920, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 142404, 141816, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 143468, 142712, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 144620, 143864, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 145961, 145016, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 147401, 146456, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 149051, 147896, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 150811, 149656, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 152802, 151416, 66, 1, nmax);

        simdtrf::transform_h_inner(buffer, 154914, 153528, 66, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 155640, 140017, 141508, 11,
                                             nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 156333, 140689, 142404, 11,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 157026, 141508, 143468, 11,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 157950, 142404, 144620, 11,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 158874, 143468, 145961, 11,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 160062, 144620, 147401, 11,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 161250, 145961, 149051, 11,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 162735, 147401, 150811, 11,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 164220, 149051, 152802, 11,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 166035, 150811, 154914, 11,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 167850, 155640, 157026, 11,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 169236, 156333, 157950, 11,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 170622, 157026, 158874, 11,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 172470, 157950, 160062, 11,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 174318, 158874, 161250, 11,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 176694, 160062, 162735, 11,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 179070, 161250, 164220, 11,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 182040, 162735, 166035, 11,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 185010, 167850, 170622, 11,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 187320, 169236, 172470, 11,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 189630, 170622, 174318, 11,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 192710, 172470, 176694, 11,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 195790, 174318, 179070, 11,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 199750, 176694, 182040, 11,
                                             nmax);

        simdtrf::compute_hrr_hg_out_of_first(buffer, coordinates, 203710, 185010, 189630, 11,
                                             nmax);

        simdtrf::compute_hrr_hg_out_of_first(buffer, coordinates, 207175, 187320, 192710, 11,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 210640, 189630, 195790, 11,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 215260, 192710, 199750, 11,
                                             nmax);

        simdtrf::compute_hrr_hh_out_of_first(buffer, coordinates, 219880, 203710, 210640, 11,
                                             nmax);

        simdtrf::compute_hrr_hh_out_of_first(buffer, coordinates, 224731, 207175, 215260, 11,
                                             nmax);

        simdtrf::transform_h_inner(buffer, 229582, 224731, 21, 11, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 229582, 121, nmax);

        simdtrf::transform_h_inner(buffer, 229582, 219880, 21, 11, nmax);

        simdtrf::transform_h_outer(values + 1331 * nvalues + n * npairs, nvalues, buffer, 229582,
                                   121, nmax);
    }

    for (size_t m = 0; m < 2662; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
