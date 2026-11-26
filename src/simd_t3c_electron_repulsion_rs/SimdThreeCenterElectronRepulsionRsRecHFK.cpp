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


#include "SimdThreeCenterElectronRepulsionRsRecHFK.hpp"

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
#include "SimdTransferHP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_hfk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_hfk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 209843, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2310 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 209843, 171608, 12585, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2978, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2981, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2984, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2987, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2990, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2993, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2996, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2999, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3002, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3005, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3008, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3011, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3014, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3017, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3020, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3023, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3026, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3029, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3032, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3035, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3038, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3041, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3044, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3047, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3050, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3053, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3056, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3059, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3062, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3065, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3068, 3, 7, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3077, 3, 8, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3086, 3, 9, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3095, 3, 10, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3104, 3, 11, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3113, 3, 12, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3122, 3, 13, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3131, 3, 14, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3140, 3, 15, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3149, 3, 16, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3158, 3, 17, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3167, 3, 18, 71,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3176, 3, 19, 74,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3185, 3, 20, 77,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3194, 3, 23, 80,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3203, 3, 24, 83,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3212, 3, 25, 86,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3221, 3, 26, 89,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3230, 3, 27, 92,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3239, 3, 28, 95,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3248, 3, 29, 98,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3257, 3, 30, 101,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3266, 3, 31, 104,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3275, 3, 32, 107,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3284, 3, 33, 110,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3293, 3, 34, 113,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3302, 3, 35, 116,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3311, 3, 36, 119,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3320, 3, 38, 122,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3338, 3, 41, 128,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3356, 3, 44, 134,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3374, 3, 47, 140,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3392, 3, 50, 146,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3410, 3, 53, 152,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3428, 3, 56, 158,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3446, 3, 59, 164,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3464, 3, 62, 170,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3482, 3, 65, 176,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3500, 3, 68, 182,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3518, 3, 71, 188,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3536, 3, 74, 194,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3554, 3, 80, 200,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3572, 3, 83, 206,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3590, 3, 86, 212,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3608, 3, 89, 218,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3626, 3, 92, 224,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3644, 3, 95, 230,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3662, 3, 98, 236,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3680, 3, 101, 242,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3698, 3, 104, 248,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3716, 3, 107, 254,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3734, 3, 110, 260,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3752, 3, 113, 266,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3770, 3, 116, 272,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3788, 3, 122, 278,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3818, 3, 128, 288,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3848, 3, 134, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3878, 3, 140, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3908, 3, 146, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3938, 3, 152, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3968, 3, 158, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3998, 3, 164, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4028, 3, 170, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4058, 3, 176, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4088, 3, 182, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4118, 3, 188, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4148, 3, 200, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4178, 3, 206, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4208, 3, 212, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4238, 3, 218, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4268, 3, 224, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4298, 3, 230, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4328, 3, 236, 458,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4358, 3, 242, 468,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4388, 3, 248, 478,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4418, 3, 254, 488,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4448, 3, 260, 498,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4478, 3, 266, 508,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4508, 3, 278, 518,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4553, 3, 288, 533,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4598, 3, 298, 548,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4643, 3, 308, 563,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4688, 3, 318, 578,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4733, 3, 328, 593,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4778, 3, 338, 608,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4823, 3, 348, 623,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4868, 3, 358, 638,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4913, 3, 368, 653,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4958, 3, 378, 668,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5003, 3, 398, 683,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5048, 3, 408, 698,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5093, 3, 418, 713,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5138, 3, 428, 728,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5183, 3, 438, 743,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5228, 3, 448, 758,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5273, 3, 458, 773,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5318, 3, 468, 788,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5363, 3, 478, 803,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5408, 3, 488, 818,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5453, 3, 498, 833,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5498, 3, 518, 848,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5561, 3, 533, 869,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5624, 3, 548, 890,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5687, 3, 563, 911,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5750, 3, 578, 932,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5813, 3, 593, 953,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5876, 3, 608, 974,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5939, 3, 623, 995,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6002, 3, 638,
                                                                       1016, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6065, 3, 653,
                                                                       1037, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6128, 3, 683,
                                                                       1058, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6191, 3, 698,
                                                                       1079, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6254, 3, 713,
                                                                       1100, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6317, 3, 728,
                                                                       1121, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6380, 3, 743,
                                                                       1142, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6443, 3, 758,
                                                                       1163, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6506, 3, 773,
                                                                       1184, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6569, 3, 788,
                                                                       1205, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6632, 3, 803,
                                                                       1226, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6695, 3, 818,
                                                                       1247, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6758, 3, 848,
                                                                       1268, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6842, 3, 869,
                                                                       1296, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6926, 3, 890,
                                                                       1324, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7010, 3, 911,
                                                                       1352, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7094, 3, 932,
                                                                       1380, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7178, 3, 953,
                                                                       1408, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7262, 3, 974,
                                                                       1436, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7346, 3, 995,
                                                                       1464, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7430, 3, 1016,
                                                                       1492, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7514, 3, 1058,
                                                                       1520, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7598, 3, 1079,
                                                                       1548, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7682, 3, 1100,
                                                                       1576, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7766, 3, 1121,
                                                                       1604, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7850, 3, 1142,
                                                                       1632, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7934, 3, 1163,
                                                                       1660, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8018, 3, 1184,
                                                                       1688, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8102, 3, 1205,
                                                                       1716, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8186, 3, 1226,
                                                                       1744, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8270, 3, 1268,
                                                                       1772, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8378, 3, 1296,
                                                                       1808, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8486, 3, 1324,
                                                                       1844, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8594, 3, 1352,
                                                                       1880, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8702, 3, 1380,
                                                                       1916, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8810, 3, 1408,
                                                                       1952, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8918, 3, 1436,
                                                                       1988, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9026, 3, 1464,
                                                                       2024, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9134, 3, 1520,
                                                                       2060, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9242, 3, 1548,
                                                                       2096, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9350, 3, 1576,
                                                                       2132, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9458, 3, 1604,
                                                                       2168, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9566, 3, 1632,
                                                                       2204, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9674, 3, 1660,
                                                                       2240, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9782, 3, 1688,
                                                                       2276, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9890, 3, 1716,
                                                                       2312, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9998, 3, 1772,
                                                                       2348, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10133, 3, 1808,
                                                                       2393, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10268, 3, 1844,
                                                                       2438, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10403, 3, 1880,
                                                                       2483, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10538, 3, 1916,
                                                                       2528, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10673, 3, 1952,
                                                                       2573, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10808, 3, 1988,
                                                                       2618, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10943, 3, 2060,
                                                                       2663, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11078, 3, 2096,
                                                                       2708, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11213, 3, 2132,
                                                                       2753, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11348, 3, 2168,
                                                                       2798, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11483, 3, 2204,
                                                                       2843, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11618, 3, 2240,
                                                                       2888, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11753, 3, 2276,
                                                                       2933, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11888, 3, 7, 8,
                                                                       2984, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11894, 3, 8, 9,
                                                                       2987, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11900, 3, 9, 10,
                                                                       2990, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11906, 3, 10, 11,
                                                                       2993, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11912, 3, 11, 12,
                                                                       2996, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11918, 3, 12, 13,
                                                                       2999, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11924, 3, 13, 14,
                                                                       3002, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11930, 3, 14, 15,
                                                                       3005, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11936, 3, 15, 16,
                                                                       3008, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11942, 3, 16, 17,
                                                                       3011, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11948, 3, 17, 18,
                                                                       3014, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11954, 3, 18, 19,
                                                                       3017, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11960, 3, 19, 20,
                                                                       3020, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11966, 3, 23, 24,
                                                                       3029, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11972, 3, 24, 25,
                                                                       3032, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11978, 3, 25, 26,
                                                                       3035, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11984, 3, 26, 27,
                                                                       3038, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11990, 3, 27, 28,
                                                                       3041, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11996, 3, 28, 29,
                                                                       3044, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12002, 3, 29, 30,
                                                                       3047, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12008, 3, 30, 31,
                                                                       3050, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12014, 3, 31, 32,
                                                                       3053, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12020, 3, 32, 33,
                                                                       3056, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12026, 3, 33, 34,
                                                                       3059, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12032, 3, 34, 35,
                                                                       3062, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12038, 3, 35, 36,
                                                                       3065, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12044, 0, 3,
                                                                       11888, 2984, 11894, 3086,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12062, 0, 3,
                                                                       11894, 2987, 11900, 3095,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12080, 0, 3,
                                                                       11900, 2990, 11906, 3104,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12098, 0, 3,
                                                                       11906, 2993, 11912, 3113,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12116, 0, 3,
                                                                       11912, 2996, 11918, 3122,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12134, 0, 3,
                                                                       11918, 2999, 11924, 3131,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12152, 0, 3,
                                                                       11924, 3002, 11930, 3140,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12170, 0, 3,
                                                                       11930, 3005, 11936, 3149,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12188, 0, 3,
                                                                       11936, 3008, 11942, 3158,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12206, 0, 3,
                                                                       11942, 3011, 11948, 3167,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12224, 0, 3,
                                                                       11948, 3014, 11954, 3176,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12242, 0, 3,
                                                                       11954, 3017, 11960, 3185,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12260, 0, 3,
                                                                       11966, 3029, 11972, 3212,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12278, 0, 3,
                                                                       11972, 3032, 11978, 3221,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12296, 0, 3,
                                                                       11978, 3035, 11984, 3230,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12314, 0, 3,
                                                                       11984, 3038, 11990, 3239,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12332, 0, 3,
                                                                       11990, 3041, 11996, 3248,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12350, 0, 3,
                                                                       11996, 3044, 12002, 3257,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12368, 0, 3,
                                                                       12002, 3047, 12008, 3266,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12386, 0, 3,
                                                                       12008, 3050, 12014, 3275,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12404, 0, 3,
                                                                       12014, 3053, 12020, 3284,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12422, 0, 3,
                                                                       12020, 3056, 12026, 3293,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12440, 0, 3,
                                                                       12026, 3059, 12032, 3302,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 12458, 0, 3,
                                                                       12032, 3062, 12038, 3311,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12476, 0, 3,
                                                                       12044, 3086, 12062, 122,
                                                                       128, 3356, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12512, 0, 3,
                                                                       12062, 3095, 12080, 128,
                                                                       134, 3374, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12548, 0, 3,
                                                                       12080, 3104, 12098, 134,
                                                                       140, 3392, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12584, 0, 3,
                                                                       12098, 3113, 12116, 140,
                                                                       146, 3410, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12620, 0, 3,
                                                                       12116, 3122, 12134, 146,
                                                                       152, 3428, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12656, 0, 3,
                                                                       12134, 3131, 12152, 152,
                                                                       158, 3446, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12692, 0, 3,
                                                                       12152, 3140, 12170, 158,
                                                                       164, 3464, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12728, 0, 3,
                                                                       12170, 3149, 12188, 164,
                                                                       170, 3482, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12764, 0, 3,
                                                                       12188, 3158, 12206, 170,
                                                                       176, 3500, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12800, 0, 3,
                                                                       12206, 3167, 12224, 176,
                                                                       182, 3518, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12836, 0, 3,
                                                                       12224, 3176, 12242, 182,
                                                                       188, 3536, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12872, 0, 3,
                                                                       12260, 3212, 12278, 200,
                                                                       206, 3590, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12908, 0, 3,
                                                                       12278, 3221, 12296, 206,
                                                                       212, 3608, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12944, 0, 3,
                                                                       12296, 3230, 12314, 212,
                                                                       218, 3626, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12980, 0, 3,
                                                                       12314, 3239, 12332, 218,
                                                                       224, 3644, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13016, 0, 3,
                                                                       12332, 3248, 12350, 224,
                                                                       230, 3662, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13052, 0, 3,
                                                                       12350, 3257, 12368, 230,
                                                                       236, 3680, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13088, 0, 3,
                                                                       12368, 3266, 12386, 236,
                                                                       242, 3698, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13124, 0, 3,
                                                                       12386, 3275, 12404, 242,
                                                                       248, 3716, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13160, 0, 3,
                                                                       12404, 3284, 12422, 248,
                                                                       254, 3734, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13196, 0, 3,
                                                                       12422, 3293, 12440, 254,
                                                                       260, 3752, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 13232, 0, 3,
                                                                       12440, 3302, 12458, 260,
                                                                       266, 3770, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13268, 0, 3,
                                                                       12476, 3356, 12512, 278,
                                                                       288, 3848, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13328, 0, 3,
                                                                       12512, 3374, 12548, 288,
                                                                       298, 3878, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13388, 0, 3,
                                                                       12548, 3392, 12584, 298,
                                                                       308, 3908, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13448, 0, 3,
                                                                       12584, 3410, 12620, 308,
                                                                       318, 3938, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13508, 0, 3,
                                                                       12620, 3428, 12656, 318,
                                                                       328, 3968, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13568, 0, 3,
                                                                       12656, 3446, 12692, 328,
                                                                       338, 3998, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13628, 0, 3,
                                                                       12692, 3464, 12728, 338,
                                                                       348, 4028, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13688, 0, 3,
                                                                       12728, 3482, 12764, 348,
                                                                       358, 4058, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13748, 0, 3,
                                                                       12764, 3500, 12800, 358,
                                                                       368, 4088, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13808, 0, 3,
                                                                       12800, 3518, 12836, 368,
                                                                       378, 4118, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13868, 0, 3,
                                                                       12872, 3590, 12908, 398,
                                                                       408, 4208, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13928, 0, 3,
                                                                       12908, 3608, 12944, 408,
                                                                       418, 4238, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13988, 0, 3,
                                                                       12944, 3626, 12980, 418,
                                                                       428, 4268, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14048, 0, 3,
                                                                       12980, 3644, 13016, 428,
                                                                       438, 4298, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14108, 0, 3,
                                                                       13016, 3662, 13052, 438,
                                                                       448, 4328, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14168, 0, 3,
                                                                       13052, 3680, 13088, 448,
                                                                       458, 4358, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14228, 0, 3,
                                                                       13088, 3698, 13124, 458,
                                                                       468, 4388, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14288, 0, 3,
                                                                       13124, 3716, 13160, 468,
                                                                       478, 4418, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14348, 0, 3,
                                                                       13160, 3734, 13196, 478,
                                                                       488, 4448, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 14408, 0, 3,
                                                                       13196, 3752, 13232, 488,
                                                                       498, 4478, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14468, 0, 3,
                                                                       13268, 3848, 13328, 518,
                                                                       533, 4598, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14558, 0, 3,
                                                                       13328, 3878, 13388, 533,
                                                                       548, 4643, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14648, 0, 3,
                                                                       13388, 3908, 13448, 548,
                                                                       563, 4688, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14738, 0, 3,
                                                                       13448, 3938, 13508, 563,
                                                                       578, 4733, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14828, 0, 3,
                                                                       13508, 3968, 13568, 578,
                                                                       593, 4778, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14918, 0, 3,
                                                                       13568, 3998, 13628, 593,
                                                                       608, 4823, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15008, 0, 3,
                                                                       13628, 4028, 13688, 608,
                                                                       623, 4868, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15098, 0, 3,
                                                                       13688, 4058, 13748, 623,
                                                                       638, 4913, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15188, 0, 3,
                                                                       13748, 4088, 13808, 638,
                                                                       653, 4958, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15278, 0, 3,
                                                                       13868, 4208, 13928, 683,
                                                                       698, 5093, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15368, 0, 3,
                                                                       13928, 4238, 13988, 698,
                                                                       713, 5138, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15458, 0, 3,
                                                                       13988, 4268, 14048, 713,
                                                                       728, 5183, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15548, 0, 3,
                                                                       14048, 4298, 14108, 728,
                                                                       743, 5228, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15638, 0, 3,
                                                                       14108, 4328, 14168, 743,
                                                                       758, 5273, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15728, 0, 3,
                                                                       14168, 4358, 14228, 758,
                                                                       773, 5318, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15818, 0, 3,
                                                                       14228, 4388, 14288, 773,
                                                                       788, 5363, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15908, 0, 3,
                                                                       14288, 4418, 14348, 788,
                                                                       803, 5408, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15998, 0, 3,
                                                                       14348, 4448, 14408, 803,
                                                                       818, 5453, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16088, 0, 3,
                                                                       14468, 4598, 14558, 848,
                                                                       869, 5624, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16214, 0, 3,
                                                                       14558, 4643, 14648, 869,
                                                                       890, 5687, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16340, 0, 3,
                                                                       14648, 4688, 14738, 890,
                                                                       911, 5750, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16466, 0, 3,
                                                                       14738, 4733, 14828, 911,
                                                                       932, 5813, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16592, 0, 3,
                                                                       14828, 4778, 14918, 932,
                                                                       953, 5876, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16718, 0, 3,
                                                                       14918, 4823, 15008, 953,
                                                                       974, 5939, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16844, 0, 3,
                                                                       15008, 4868, 15098, 974,
                                                                       995, 6002, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 16970, 0, 3,
                                                                       15098, 4913, 15188, 995,
                                                                       1016, 6065, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17096, 0, 3,
                                                                       15278, 5093, 15368, 1058,
                                                                       1079, 6254, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17222, 0, 3,
                                                                       15368, 5138, 15458, 1079,
                                                                       1100, 6317, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17348, 0, 3,
                                                                       15458, 5183, 15548, 1100,
                                                                       1121, 6380, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17474, 0, 3,
                                                                       15548, 5228, 15638, 1121,
                                                                       1142, 6443, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17600, 0, 3,
                                                                       15638, 5273, 15728, 1142,
                                                                       1163, 6506, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17726, 0, 3,
                                                                       15728, 5318, 15818, 1163,
                                                                       1184, 6569, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17852, 0, 3,
                                                                       15818, 5363, 15908, 1184,
                                                                       1205, 6632, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17978, 0, 3,
                                                                       15908, 5408, 15998, 1205,
                                                                       1226, 6695, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18104, 0, 3,
                                                                       16088, 5624, 16214, 1268,
                                                                       1296, 6926, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18272, 0, 3,
                                                                       16214, 5687, 16340, 1296,
                                                                       1324, 7010, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18440, 0, 3,
                                                                       16340, 5750, 16466, 1324,
                                                                       1352, 7094, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18608, 0, 3,
                                                                       16466, 5813, 16592, 1352,
                                                                       1380, 7178, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18776, 0, 3,
                                                                       16592, 5876, 16718, 1380,
                                                                       1408, 7262, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18944, 0, 3,
                                                                       16718, 5939, 16844, 1408,
                                                                       1436, 7346, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19112, 0, 3,
                                                                       16844, 6002, 16970, 1436,
                                                                       1464, 7430, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19280, 0, 3,
                                                                       17096, 6254, 17222, 1520,
                                                                       1548, 7682, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19448, 0, 3,
                                                                       17222, 6317, 17348, 1548,
                                                                       1576, 7766, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19616, 0, 3,
                                                                       17348, 6380, 17474, 1576,
                                                                       1604, 7850, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19784, 0, 3,
                                                                       17474, 6443, 17600, 1604,
                                                                       1632, 7934, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19952, 0, 3,
                                                                       17600, 6506, 17726, 1632,
                                                                       1660, 8018, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20120, 0, 3,
                                                                       17726, 6569, 17852, 1660,
                                                                       1688, 8102, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20288, 0, 3,
                                                                       17852, 6632, 17978, 1688,
                                                                       1716, 8186, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20456, 0, 3,
                                                                       18104, 6926, 18272, 1772,
                                                                       1808, 8486, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20672, 0, 3,
                                                                       18272, 7010, 18440, 1808,
                                                                       1844, 8594, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20888, 0, 3,
                                                                       18440, 7094, 18608, 1844,
                                                                       1880, 8702, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21104, 0, 3,
                                                                       18608, 7178, 18776, 1880,
                                                                       1916, 8810, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21320, 0, 3,
                                                                       18776, 7262, 18944, 1916,
                                                                       1952, 8918, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21536, 0, 3,
                                                                       18944, 7346, 19112, 1952,
                                                                       1988, 9026, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21752, 0, 3,
                                                                       19280, 7682, 19448, 2060,
                                                                       2096, 9350, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21968, 0, 3,
                                                                       19448, 7766, 19616, 2096,
                                                                       2132, 9458, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 22184, 0, 3,
                                                                       19616, 7850, 19784, 2132,
                                                                       2168, 9566, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 22400, 0, 3,
                                                                       19784, 7934, 19952, 2168,
                                                                       2204, 9674, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 22616, 0, 3,
                                                                       19952, 8018, 20120, 2204,
                                                                       2240, 9782, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 22832, 0, 3,
                                                                       20120, 8102, 20288, 2240,
                                                                       2276, 9890, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23048, 0, 3,
                                                                       20456, 8486, 20672, 2348,
                                                                       2393, 10268, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23318, 0, 3,
                                                                       20672, 8594, 20888, 2393,
                                                                       2438, 10403, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23588, 0, 3,
                                                                       20888, 8702, 21104, 2438,
                                                                       2483, 10538, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23858, 0, 3,
                                                                       21104, 8810, 21320, 2483,
                                                                       2528, 10673, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 24128, 0, 3,
                                                                       21320, 8918, 21536, 2528,
                                                                       2573, 10808, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 24398, 0, 3,
                                                                       21752, 9350, 21968, 2663,
                                                                       2708, 11213, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 24668, 0, 3,
                                                                       21968, 9458, 22184, 2708,
                                                                       2753, 11348, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 24938, 0, 3,
                                                                       22184, 9566, 22400, 2753,
                                                                       2798, 11483, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 25208, 0, 3,
                                                                       22400, 9674, 22616, 2798,
                                                                       2843, 11618, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 25478, 0, 3,
                                                                       22616, 9782, 22832, 2843,
                                                                       2888, 11753, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25748, 3, 2978,
                                                                       2981, 11888, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25758, 3, 2981,
                                                                       2984, 11894, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25768, 3, 2984,
                                                                       2987, 11900, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25778, 3, 2987,
                                                                       2990, 11906, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25788, 3, 2990,
                                                                       2993, 11912, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25798, 3, 2993,
                                                                       2996, 11918, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25808, 3, 2996,
                                                                       2999, 11924, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25818, 3, 2999,
                                                                       3002, 11930, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25828, 3, 3002,
                                                                       3005, 11936, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25838, 3, 3005,
                                                                       3008, 11942, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25848, 3, 3008,
                                                                       3011, 11948, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25858, 3, 3011,
                                                                       3014, 11954, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25868, 3, 3014,
                                                                       3017, 11960, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25878, 3, 3023,
                                                                       3026, 11966, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25888, 3, 3026,
                                                                       3029, 11972, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25898, 3, 3029,
                                                                       3032, 11978, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25908, 3, 3032,
                                                                       3035, 11984, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25918, 3, 3035,
                                                                       3038, 11990, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25928, 3, 3038,
                                                                       3041, 11996, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25938, 3, 3041,
                                                                       3044, 12002, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25948, 3, 3044,
                                                                       3047, 12008, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25958, 3, 3047,
                                                                       3050, 12014, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25968, 3, 3050,
                                                                       3053, 12020, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25978, 3, 3053,
                                                                       3056, 12026, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25988, 3, 3056,
                                                                       3059, 12032, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 25998, 3, 3059,
                                                                       3062, 12038, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26008, 0, 3,
                                                                       25748, 11888, 25758, 3068,
                                                                       3077, 12044, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26038, 0, 3,
                                                                       25758, 11894, 25768, 3077,
                                                                       3086, 12062, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26068, 0, 3,
                                                                       25768, 11900, 25778, 3086,
                                                                       3095, 12080, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26098, 0, 3,
                                                                       25778, 11906, 25788, 3095,
                                                                       3104, 12098, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26128, 0, 3,
                                                                       25788, 11912, 25798, 3104,
                                                                       3113, 12116, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26158, 0, 3,
                                                                       25798, 11918, 25808, 3113,
                                                                       3122, 12134, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26188, 0, 3,
                                                                       25808, 11924, 25818, 3122,
                                                                       3131, 12152, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26218, 0, 3,
                                                                       25818, 11930, 25828, 3131,
                                                                       3140, 12170, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26248, 0, 3,
                                                                       25828, 11936, 25838, 3140,
                                                                       3149, 12188, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26278, 0, 3,
                                                                       25838, 11942, 25848, 3149,
                                                                       3158, 12206, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26308, 0, 3,
                                                                       25848, 11948, 25858, 3158,
                                                                       3167, 12224, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26338, 0, 3,
                                                                       25858, 11954, 25868, 3167,
                                                                       3176, 12242, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26368, 0, 3,
                                                                       25878, 11966, 25888, 3194,
                                                                       3203, 12260, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26398, 0, 3,
                                                                       25888, 11972, 25898, 3203,
                                                                       3212, 12278, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26428, 0, 3,
                                                                       25898, 11978, 25908, 3212,
                                                                       3221, 12296, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26458, 0, 3,
                                                                       25908, 11984, 25918, 3221,
                                                                       3230, 12314, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26488, 0, 3,
                                                                       25918, 11990, 25928, 3230,
                                                                       3239, 12332, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26518, 0, 3,
                                                                       25928, 11996, 25938, 3239,
                                                                       3248, 12350, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26548, 0, 3,
                                                                       25938, 12002, 25948, 3248,
                                                                       3257, 12368, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26578, 0, 3,
                                                                       25948, 12008, 25958, 3257,
                                                                       3266, 12386, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26608, 0, 3,
                                                                       25958, 12014, 25968, 3266,
                                                                       3275, 12404, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26638, 0, 3,
                                                                       25968, 12020, 25978, 3275,
                                                                       3284, 12422, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26668, 0, 3,
                                                                       25978, 12026, 25988, 3284,
                                                                       3293, 12440, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 26698, 0, 3,
                                                                       25988, 12032, 25998, 3293,
                                                                       3302, 12458, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 26728, 0, 3,
                                                                       26008, 12044, 26038, 3320,
                                                                       3338, 12476, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 26788, 0, 3,
                                                                       26038, 12062, 26068, 3338,
                                                                       3356, 12512, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 26848, 0, 3,
                                                                       26068, 12080, 26098, 3356,
                                                                       3374, 12548, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 26908, 0, 3,
                                                                       26098, 12098, 26128, 3374,
                                                                       3392, 12584, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 26968, 0, 3,
                                                                       26128, 12116, 26158, 3392,
                                                                       3410, 12620, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27028, 0, 3,
                                                                       26158, 12134, 26188, 3410,
                                                                       3428, 12656, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27088, 0, 3,
                                                                       26188, 12152, 26218, 3428,
                                                                       3446, 12692, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27148, 0, 3,
                                                                       26218, 12170, 26248, 3446,
                                                                       3464, 12728, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27208, 0, 3,
                                                                       26248, 12188, 26278, 3464,
                                                                       3482, 12764, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27268, 0, 3,
                                                                       26278, 12206, 26308, 3482,
                                                                       3500, 12800, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27328, 0, 3,
                                                                       26308, 12224, 26338, 3500,
                                                                       3518, 12836, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27388, 0, 3,
                                                                       26368, 12260, 26398, 3554,
                                                                       3572, 12872, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27448, 0, 3,
                                                                       26398, 12278, 26428, 3572,
                                                                       3590, 12908, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27508, 0, 3,
                                                                       26428, 12296, 26458, 3590,
                                                                       3608, 12944, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27568, 0, 3,
                                                                       26458, 12314, 26488, 3608,
                                                                       3626, 12980, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27628, 0, 3,
                                                                       26488, 12332, 26518, 3626,
                                                                       3644, 13016, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27688, 0, 3,
                                                                       26518, 12350, 26548, 3644,
                                                                       3662, 13052, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27748, 0, 3,
                                                                       26548, 12368, 26578, 3662,
                                                                       3680, 13088, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27808, 0, 3,
                                                                       26578, 12386, 26608, 3680,
                                                                       3698, 13124, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27868, 0, 3,
                                                                       26608, 12404, 26638, 3698,
                                                                       3716, 13160, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27928, 0, 3,
                                                                       26638, 12422, 26668, 3716,
                                                                       3734, 13196, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 27988, 0, 3,
                                                                       26668, 12440, 26698, 3734,
                                                                       3752, 13232, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 28048, 0, 3,
                                                                       26728, 12476, 26788, 3788,
                                                                       3818, 13268, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 28148, 0, 3,
                                                                       26788, 12512, 26848, 3818,
                                                                       3848, 13328, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 28248, 0, 3,
                                                                       26848, 12548, 26908, 3848,
                                                                       3878, 13388, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 28348, 0, 3,
                                                                       26908, 12584, 26968, 3878,
                                                                       3908, 13448, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 28448, 0, 3,
                                                                       26968, 12620, 27028, 3908,
                                                                       3938, 13508, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 28548, 0, 3,
                                                                       27028, 12656, 27088, 3938,
                                                                       3968, 13568, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 28648, 0, 3,
                                                                       27088, 12692, 27148, 3968,
                                                                       3998, 13628, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 28748, 0, 3,
                                                                       27148, 12728, 27208, 3998,
                                                                       4028, 13688, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 28848, 0, 3,
                                                                       27208, 12764, 27268, 4028,
                                                                       4058, 13748, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 28948, 0, 3,
                                                                       27268, 12800, 27328, 4058,
                                                                       4088, 13808, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29048, 0, 3,
                                                                       27388, 12872, 27448, 4148,
                                                                       4178, 13868, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29148, 0, 3,
                                                                       27448, 12908, 27508, 4178,
                                                                       4208, 13928, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29248, 0, 3,
                                                                       27508, 12944, 27568, 4208,
                                                                       4238, 13988, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29348, 0, 3,
                                                                       27568, 12980, 27628, 4238,
                                                                       4268, 14048, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29448, 0, 3,
                                                                       27628, 13016, 27688, 4268,
                                                                       4298, 14108, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29548, 0, 3,
                                                                       27688, 13052, 27748, 4298,
                                                                       4328, 14168, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29648, 0, 3,
                                                                       27748, 13088, 27808, 4328,
                                                                       4358, 14228, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29748, 0, 3,
                                                                       27808, 13124, 27868, 4358,
                                                                       4388, 14288, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29848, 0, 3,
                                                                       27868, 13160, 27928, 4388,
                                                                       4418, 14348, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29948, 0, 3,
                                                                       27928, 13196, 27988, 4418,
                                                                       4448, 14408, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30048, 0, 3,
                                                                       28048, 13268, 28148, 4508,
                                                                       4553, 14468, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30198, 0, 3,
                                                                       28148, 13328, 28248, 4553,
                                                                       4598, 14558, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30348, 0, 3,
                                                                       28248, 13388, 28348, 4598,
                                                                       4643, 14648, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30498, 0, 3,
                                                                       28348, 13448, 28448, 4643,
                                                                       4688, 14738, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30648, 0, 3,
                                                                       28448, 13508, 28548, 4688,
                                                                       4733, 14828, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30798, 0, 3,
                                                                       28548, 13568, 28648, 4733,
                                                                       4778, 14918, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30948, 0, 3,
                                                                       28648, 13628, 28748, 4778,
                                                                       4823, 15008, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31098, 0, 3,
                                                                       28748, 13688, 28848, 4823,
                                                                       4868, 15098, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31248, 0, 3,
                                                                       28848, 13748, 28948, 4868,
                                                                       4913, 15188, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31398, 0, 3,
                                                                       29048, 13868, 29148, 5003,
                                                                       5048, 15278, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31548, 0, 3,
                                                                       29148, 13928, 29248, 5048,
                                                                       5093, 15368, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31698, 0, 3,
                                                                       29248, 13988, 29348, 5093,
                                                                       5138, 15458, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31848, 0, 3,
                                                                       29348, 14048, 29448, 5138,
                                                                       5183, 15548, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31998, 0, 3,
                                                                       29448, 14108, 29548, 5183,
                                                                       5228, 15638, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 32148, 0, 3,
                                                                       29548, 14168, 29648, 5228,
                                                                       5273, 15728, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 32298, 0, 3,
                                                                       29648, 14228, 29748, 5273,
                                                                       5318, 15818, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 32448, 0, 3,
                                                                       29748, 14288, 29848, 5318,
                                                                       5363, 15908, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 32598, 0, 3,
                                                                       29848, 14348, 29948, 5363,
                                                                       5408, 15998, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32748, 0, 3,
                                                                       30048, 14468, 30198, 5498,
                                                                       5561, 16088, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32958, 0, 3,
                                                                       30198, 14558, 30348, 5561,
                                                                       5624, 16214, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33168, 0, 3,
                                                                       30348, 14648, 30498, 5624,
                                                                       5687, 16340, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33378, 0, 3,
                                                                       30498, 14738, 30648, 5687,
                                                                       5750, 16466, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33588, 0, 3,
                                                                       30648, 14828, 30798, 5750,
                                                                       5813, 16592, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33798, 0, 3,
                                                                       30798, 14918, 30948, 5813,
                                                                       5876, 16718, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 34008, 0, 3,
                                                                       30948, 15008, 31098, 5876,
                                                                       5939, 16844, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 34218, 0, 3,
                                                                       31098, 15098, 31248, 5939,
                                                                       6002, 16970, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 34428, 0, 3,
                                                                       31398, 15278, 31548, 6128,
                                                                       6191, 17096, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 34638, 0, 3,
                                                                       31548, 15368, 31698, 6191,
                                                                       6254, 17222, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 34848, 0, 3,
                                                                       31698, 15458, 31848, 6254,
                                                                       6317, 17348, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35058, 0, 3,
                                                                       31848, 15548, 31998, 6317,
                                                                       6380, 17474, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35268, 0, 3,
                                                                       31998, 15638, 32148, 6380,
                                                                       6443, 17600, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35478, 0, 3,
                                                                       32148, 15728, 32298, 6443,
                                                                       6506, 17726, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35688, 0, 3,
                                                                       32298, 15818, 32448, 6506,
                                                                       6569, 17852, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 35898, 0, 3,
                                                                       32448, 15908, 32598, 6569,
                                                                       6632, 17978, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 36108, 0, 3,
                                                                       32748, 16088, 32958, 6758,
                                                                       6842, 18104, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 36388, 0, 3,
                                                                       32958, 16214, 33168, 6842,
                                                                       6926, 18272, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 36668, 0, 3,
                                                                       33168, 16340, 33378, 6926,
                                                                       7010, 18440, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 36948, 0, 3,
                                                                       33378, 16466, 33588, 7010,
                                                                       7094, 18608, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 37228, 0, 3,
                                                                       33588, 16592, 33798, 7094,
                                                                       7178, 18776, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 37508, 0, 3,
                                                                       33798, 16718, 34008, 7178,
                                                                       7262, 18944, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 37788, 0, 3,
                                                                       34008, 16844, 34218, 7262,
                                                                       7346, 19112, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38068, 0, 3,
                                                                       34428, 17096, 34638, 7514,
                                                                       7598, 19280, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38348, 0, 3,
                                                                       34638, 17222, 34848, 7598,
                                                                       7682, 19448, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38628, 0, 3,
                                                                       34848, 17348, 35058, 7682,
                                                                       7766, 19616, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 38908, 0, 3,
                                                                       35058, 17474, 35268, 7766,
                                                                       7850, 19784, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 39188, 0, 3,
                                                                       35268, 17600, 35478, 7850,
                                                                       7934, 19952, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 39468, 0, 3,
                                                                       35478, 17726, 35688, 7934,
                                                                       8018, 20120, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 39748, 0, 3,
                                                                       35688, 17852, 35898, 8018,
                                                                       8102, 20288, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 40028, 0, 3,
                                                                       36108, 18104, 36388, 8270,
                                                                       8378, 20456, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 40388, 0, 3,
                                                                       36388, 18272, 36668, 8378,
                                                                       8486, 20672, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 40748, 0, 3,
                                                                       36668, 18440, 36948, 8486,
                                                                       8594, 20888, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 41108, 0, 3,
                                                                       36948, 18608, 37228, 8594,
                                                                       8702, 21104, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 41468, 0, 3,
                                                                       37228, 18776, 37508, 8702,
                                                                       8810, 21320, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 41828, 0, 3,
                                                                       37508, 18944, 37788, 8810,
                                                                       8918, 21536, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 42188, 0, 3,
                                                                       38068, 19280, 38348, 9134,
                                                                       9242, 21752, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 42548, 0, 3,
                                                                       38348, 19448, 38628, 9242,
                                                                       9350, 21968, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 42908, 0, 3,
                                                                       38628, 19616, 38908, 9350,
                                                                       9458, 22184, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 43268, 0, 3,
                                                                       38908, 19784, 39188, 9458,
                                                                       9566, 22400, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 43628, 0, 3,
                                                                       39188, 19952, 39468, 9566,
                                                                       9674, 22616, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 43988, 0, 3,
                                                                       39468, 20120, 39748, 9674,
                                                                       9782, 22832, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 44348, 0, 3,
                                                                       40028, 20456, 40388, 9998,
                                                                       10133, 23048, ncols,
                                                                       gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 44798, 0, 3,
                                                                       40388, 20672, 40748,
                                                                       10133, 10268, 23318,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 45248, 0, 3,
                                                                       40748, 20888, 41108,
                                                                       10268, 10403, 23588,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 45698, 0, 3,
                                                                       41108, 21104, 41468,
                                                                       10403, 10538, 23858,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 46148, 0, 3,
                                                                       41468, 21320, 41828,
                                                                       10538, 10673, 24128,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 46598, 0, 3,
                                                                       42188, 21752, 42548,
                                                                       10943, 11078, 24398,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 47048, 0, 3,
                                                                       42548, 21968, 42908,
                                                                       11078, 11213, 24668,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 47498, 0, 3,
                                                                       42908, 22184, 43268,
                                                                       11213, 11348, 24938,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 47948, 0, 3,
                                                                       43268, 22400, 43628,
                                                                       11348, 11483, 25208,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 48398, 0, 3,
                                                                       43628, 22616, 43988,
                                                                       11483, 11618, 25478,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48848, 3, 11888,
                                                                       11894, 25768, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48863, 3, 11894,
                                                                       11900, 25778, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48878, 3, 11900,
                                                                       11906, 25788, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48893, 3, 11906,
                                                                       11912, 25798, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48908, 3, 11912,
                                                                       11918, 25808, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48923, 3, 11918,
                                                                       11924, 25818, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48938, 3, 11924,
                                                                       11930, 25828, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48953, 3, 11930,
                                                                       11936, 25838, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48968, 3, 11936,
                                                                       11942, 25848, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48983, 3, 11942,
                                                                       11948, 25858, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 48998, 3, 11948,
                                                                       11954, 25868, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49013, 3, 11966,
                                                                       11972, 25898, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49028, 3, 11972,
                                                                       11978, 25908, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49043, 3, 11978,
                                                                       11984, 25918, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49058, 3, 11984,
                                                                       11990, 25928, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49073, 3, 11990,
                                                                       11996, 25938, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49088, 3, 11996,
                                                                       12002, 25948, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49103, 3, 12002,
                                                                       12008, 25958, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49118, 3, 12008,
                                                                       12014, 25968, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49133, 3, 12014,
                                                                       12020, 25978, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49148, 3, 12020,
                                                                       12026, 25988, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49163, 3, 12026,
                                                                       12032, 25998, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49178, 0, 3,
                                                                       48848, 25768, 48863,
                                                                       12044, 12062, 26068,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49223, 0, 3,
                                                                       48863, 25778, 48878,
                                                                       12062, 12080, 26098,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49268, 0, 3,
                                                                       48878, 25788, 48893,
                                                                       12080, 12098, 26128,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49313, 0, 3,
                                                                       48893, 25798, 48908,
                                                                       12098, 12116, 26158,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49358, 0, 3,
                                                                       48908, 25808, 48923,
                                                                       12116, 12134, 26188,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49403, 0, 3,
                                                                       48923, 25818, 48938,
                                                                       12134, 12152, 26218,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49448, 0, 3,
                                                                       48938, 25828, 48953,
                                                                       12152, 12170, 26248,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49493, 0, 3,
                                                                       48953, 25838, 48968,
                                                                       12170, 12188, 26278,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49538, 0, 3,
                                                                       48968, 25848, 48983,
                                                                       12188, 12206, 26308,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49583, 0, 3,
                                                                       48983, 25858, 48998,
                                                                       12206, 12224, 26338,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49628, 0, 3,
                                                                       49013, 25898, 49028,
                                                                       12260, 12278, 26428,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49673, 0, 3,
                                                                       49028, 25908, 49043,
                                                                       12278, 12296, 26458,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49718, 0, 3,
                                                                       49043, 25918, 49058,
                                                                       12296, 12314, 26488,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49763, 0, 3,
                                                                       49058, 25928, 49073,
                                                                       12314, 12332, 26518,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49808, 0, 3,
                                                                       49073, 25938, 49088,
                                                                       12332, 12350, 26548,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49853, 0, 3,
                                                                       49088, 25948, 49103,
                                                                       12350, 12368, 26578,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49898, 0, 3,
                                                                       49103, 25958, 49118,
                                                                       12368, 12386, 26608,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49943, 0, 3,
                                                                       49118, 25968, 49133,
                                                                       12386, 12404, 26638,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49988, 0, 3,
                                                                       49133, 25978, 49148,
                                                                       12404, 12422, 26668,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50033, 0, 3,
                                                                       49148, 25988, 49163,
                                                                       12422, 12440, 26698,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50078, 0, 3,
                                                                       49178, 26068, 49223,
                                                                       12476, 12512, 26848,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50168, 0, 3,
                                                                       49223, 26098, 49268,
                                                                       12512, 12548, 26908,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50258, 0, 3,
                                                                       49268, 26128, 49313,
                                                                       12548, 12584, 26968,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50348, 0, 3,
                                                                       49313, 26158, 49358,
                                                                       12584, 12620, 27028,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50438, 0, 3,
                                                                       49358, 26188, 49403,
                                                                       12620, 12656, 27088,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50528, 0, 3,
                                                                       49403, 26218, 49448,
                                                                       12656, 12692, 27148,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50618, 0, 3,
                                                                       49448, 26248, 49493,
                                                                       12692, 12728, 27208,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50708, 0, 3,
                                                                       49493, 26278, 49538,
                                                                       12728, 12764, 27268,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50798, 0, 3,
                                                                       49538, 26308, 49583,
                                                                       12764, 12800, 27328,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50888, 0, 3,
                                                                       49628, 26428, 49673,
                                                                       12872, 12908, 27508,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50978, 0, 3,
                                                                       49673, 26458, 49718,
                                                                       12908, 12944, 27568,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51068, 0, 3,
                                                                       49718, 26488, 49763,
                                                                       12944, 12980, 27628,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51158, 0, 3,
                                                                       49763, 26518, 49808,
                                                                       12980, 13016, 27688,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51248, 0, 3,
                                                                       49808, 26548, 49853,
                                                                       13016, 13052, 27748,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51338, 0, 3,
                                                                       49853, 26578, 49898,
                                                                       13052, 13088, 27808,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51428, 0, 3,
                                                                       49898, 26608, 49943,
                                                                       13088, 13124, 27868,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51518, 0, 3,
                                                                       49943, 26638, 49988,
                                                                       13124, 13160, 27928,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51608, 0, 3,
                                                                       49988, 26668, 50033,
                                                                       13160, 13196, 27988,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 51698, 0, 3,
                                                                       50078, 26848, 50168,
                                                                       13268, 13328, 28248,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 51848, 0, 3,
                                                                       50168, 26908, 50258,
                                                                       13328, 13388, 28348,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 51998, 0, 3,
                                                                       50258, 26968, 50348,
                                                                       13388, 13448, 28448,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52148, 0, 3,
                                                                       50348, 27028, 50438,
                                                                       13448, 13508, 28548,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52298, 0, 3,
                                                                       50438, 27088, 50528,
                                                                       13508, 13568, 28648,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52448, 0, 3,
                                                                       50528, 27148, 50618,
                                                                       13568, 13628, 28748,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52598, 0, 3,
                                                                       50618, 27208, 50708,
                                                                       13628, 13688, 28848,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52748, 0, 3,
                                                                       50708, 27268, 50798,
                                                                       13688, 13748, 28948,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52898, 0, 3,
                                                                       50888, 27508, 50978,
                                                                       13868, 13928, 29248,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 53048, 0, 3,
                                                                       50978, 27568, 51068,
                                                                       13928, 13988, 29348,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 53198, 0, 3,
                                                                       51068, 27628, 51158,
                                                                       13988, 14048, 29448,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 53348, 0, 3,
                                                                       51158, 27688, 51248,
                                                                       14048, 14108, 29548,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 53498, 0, 3,
                                                                       51248, 27748, 51338,
                                                                       14108, 14168, 29648,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 53648, 0, 3,
                                                                       51338, 27808, 51428,
                                                                       14168, 14228, 29748,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 53798, 0, 3,
                                                                       51428, 27868, 51518,
                                                                       14228, 14288, 29848,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 53948, 0, 3,
                                                                       51518, 27928, 51608,
                                                                       14288, 14348, 29948,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54098, 0, 3,
                                                                       51698, 28248, 51848,
                                                                       14468, 14558, 30348,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54323, 0, 3,
                                                                       51848, 28348, 51998,
                                                                       14558, 14648, 30498,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54548, 0, 3,
                                                                       51998, 28448, 52148,
                                                                       14648, 14738, 30648,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54773, 0, 3,
                                                                       52148, 28548, 52298,
                                                                       14738, 14828, 30798,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54998, 0, 3,
                                                                       52298, 28648, 52448,
                                                                       14828, 14918, 30948,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 55223, 0, 3,
                                                                       52448, 28748, 52598,
                                                                       14918, 15008, 31098,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 55448, 0, 3,
                                                                       52598, 28848, 52748,
                                                                       15008, 15098, 31248,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 55673, 0, 3,
                                                                       52898, 29248, 53048,
                                                                       15278, 15368, 31698,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 55898, 0, 3,
                                                                       53048, 29348, 53198,
                                                                       15368, 15458, 31848,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 56123, 0, 3,
                                                                       53198, 29448, 53348,
                                                                       15458, 15548, 31998,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 56348, 0, 3,
                                                                       53348, 29548, 53498,
                                                                       15548, 15638, 32148,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 56573, 0, 3,
                                                                       53498, 29648, 53648,
                                                                       15638, 15728, 32298,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 56798, 0, 3,
                                                                       53648, 29748, 53798,
                                                                       15728, 15818, 32448,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 57023, 0, 3,
                                                                       53798, 29848, 53948,
                                                                       15818, 15908, 32598,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 57248, 0, 3,
                                                                       54098, 30348, 54323,
                                                                       16088, 16214, 33168,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 57563, 0, 3,
                                                                       54323, 30498, 54548,
                                                                       16214, 16340, 33378,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 57878, 0, 3,
                                                                       54548, 30648, 54773,
                                                                       16340, 16466, 33588,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 58193, 0, 3,
                                                                       54773, 30798, 54998,
                                                                       16466, 16592, 33798,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 58508, 0, 3,
                                                                       54998, 30948, 55223,
                                                                       16592, 16718, 34008,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 58823, 0, 3,
                                                                       55223, 31098, 55448,
                                                                       16718, 16844, 34218,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 59138, 0, 3,
                                                                       55673, 31698, 55898,
                                                                       17096, 17222, 34848,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 59453, 0, 3,
                                                                       55898, 31848, 56123,
                                                                       17222, 17348, 35058,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 59768, 0, 3,
                                                                       56123, 31998, 56348,
                                                                       17348, 17474, 35268,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 60083, 0, 3,
                                                                       56348, 32148, 56573,
                                                                       17474, 17600, 35478,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 60398, 0, 3,
                                                                       56573, 32298, 56798,
                                                                       17600, 17726, 35688,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 60713, 0, 3,
                                                                       56798, 32448, 57023,
                                                                       17726, 17852, 35898,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 61028, 0, 3,
                                                                       57248, 33168, 57563,
                                                                       18104, 18272, 36668,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 61448, 0, 3,
                                                                       57563, 33378, 57878,
                                                                       18272, 18440, 36948,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 61868, 0, 3,
                                                                       57878, 33588, 58193,
                                                                       18440, 18608, 37228,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 62288, 0, 3,
                                                                       58193, 33798, 58508,
                                                                       18608, 18776, 37508,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 62708, 0, 3,
                                                                       58508, 34008, 58823,
                                                                       18776, 18944, 37788,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 63128, 0, 3,
                                                                       59138, 34848, 59453,
                                                                       19280, 19448, 38628,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 63548, 0, 3,
                                                                       59453, 35058, 59768,
                                                                       19448, 19616, 38908,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 63968, 0, 3,
                                                                       59768, 35268, 60083,
                                                                       19616, 19784, 39188,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 64388, 0, 3,
                                                                       60083, 35478, 60398,
                                                                       19784, 19952, 39468,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 64808, 0, 3,
                                                                       60398, 35688, 60713,
                                                                       19952, 20120, 39748,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 65228, 0, 3,
                                                                       61028, 36668, 61448,
                                                                       20456, 20672, 40748,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 65768, 0, 3,
                                                                       61448, 36948, 61868,
                                                                       20672, 20888, 41108,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 66308, 0, 3,
                                                                       61868, 37228, 62288,
                                                                       20888, 21104, 41468,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 66848, 0, 3,
                                                                       62288, 37508, 62708,
                                                                       21104, 21320, 41828,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 67388, 0, 3,
                                                                       63128, 38628, 63548,
                                                                       21752, 21968, 42908,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 67928, 0, 3,
                                                                       63548, 38908, 63968,
                                                                       21968, 22184, 43268,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 68468, 0, 3,
                                                                       63968, 39188, 64388,
                                                                       22184, 22400, 43628,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 69008, 0, 3,
                                                                       64388, 39468, 64808,
                                                                       22400, 22616, 43988,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 69548, 0, 3,
                                                                       65228, 40748, 65768,
                                                                       23048, 23318, 45248,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 70223, 0, 3,
                                                                       65768, 41108, 66308,
                                                                       23318, 23588, 45698,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 70898, 0, 3,
                                                                       66308, 41468, 66848,
                                                                       23588, 23858, 46148,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 71573, 0, 3,
                                                                       67388, 42908, 67928,
                                                                       24398, 24668, 47498,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 72248, 0, 3,
                                                                       67928, 43268, 68468,
                                                                       24668, 24938, 47948,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 72923, 0, 3,
                                                                       68468, 43628, 69008,
                                                                       24938, 25208, 48398,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73598, 3, 25748,
                                                                       25758, 48848, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73619, 3, 25758,
                                                                       25768, 48863, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73640, 3, 25768,
                                                                       25778, 48878, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73661, 3, 25778,
                                                                       25788, 48893, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73682, 3, 25788,
                                                                       25798, 48908, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73703, 3, 25798,
                                                                       25808, 48923, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73724, 3, 25808,
                                                                       25818, 48938, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73745, 3, 25818,
                                                                       25828, 48953, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73766, 3, 25828,
                                                                       25838, 48968, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73787, 3, 25838,
                                                                       25848, 48983, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73808, 3, 25848,
                                                                       25858, 48998, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73829, 3, 25878,
                                                                       25888, 49013, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73850, 3, 25888,
                                                                       25898, 49028, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73871, 3, 25898,
                                                                       25908, 49043, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73892, 3, 25908,
                                                                       25918, 49058, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73913, 3, 25918,
                                                                       25928, 49073, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73934, 3, 25928,
                                                                       25938, 49088, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73955, 3, 25938,
                                                                       25948, 49103, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73976, 3, 25948,
                                                                       25958, 49118, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 73997, 3, 25958,
                                                                       25968, 49133, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 74018, 3, 25968,
                                                                       25978, 49148, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 74039, 3, 25978,
                                                                       25988, 49163, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 74060, 0, 3,
                                                                       73598, 48848, 73619,
                                                                       26008, 26038, 49178,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 74123, 0, 3,
                                                                       73619, 48863, 73640,
                                                                       26038, 26068, 49223,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 74186, 0, 3,
                                                                       73640, 48878, 73661,
                                                                       26068, 26098, 49268,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 74249, 0, 3,
                                                                       73661, 48893, 73682,
                                                                       26098, 26128, 49313,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 74312, 0, 3,
                                                                       73682, 48908, 73703,
                                                                       26128, 26158, 49358,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 74375, 0, 3,
                                                                       73703, 48923, 73724,
                                                                       26158, 26188, 49403,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 74438, 0, 3,
                                                                       73724, 48938, 73745,
                                                                       26188, 26218, 49448,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 74501, 0, 3,
                                                                       73745, 48953, 73766,
                                                                       26218, 26248, 49493,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 74564, 0, 3,
                                                                       73766, 48968, 73787,
                                                                       26248, 26278, 49538,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 74627, 0, 3,
                                                                       73787, 48983, 73808,
                                                                       26278, 26308, 49583,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 74690, 0, 3,
                                                                       73829, 49013, 73850,
                                                                       26368, 26398, 49628,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 74753, 0, 3,
                                                                       73850, 49028, 73871,
                                                                       26398, 26428, 49673,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 74816, 0, 3,
                                                                       73871, 49043, 73892,
                                                                       26428, 26458, 49718,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 74879, 0, 3,
                                                                       73892, 49058, 73913,
                                                                       26458, 26488, 49763,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 74942, 0, 3,
                                                                       73913, 49073, 73934,
                                                                       26488, 26518, 49808,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 75005, 0, 3,
                                                                       73934, 49088, 73955,
                                                                       26518, 26548, 49853,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 75068, 0, 3,
                                                                       73955, 49103, 73976,
                                                                       26548, 26578, 49898,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 75131, 0, 3,
                                                                       73976, 49118, 73997,
                                                                       26578, 26608, 49943,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 75194, 0, 3,
                                                                       73997, 49133, 74018,
                                                                       26608, 26638, 49988,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 75257, 0, 3,
                                                                       74018, 49148, 74039,
                                                                       26638, 26668, 50033,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 75320, 0, 3,
                                                                       74060, 49178, 74123,
                                                                       26728, 26788, 50078,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 75446, 0, 3,
                                                                       74123, 49223, 74186,
                                                                       26788, 26848, 50168,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 75572, 0, 3,
                                                                       74186, 49268, 74249,
                                                                       26848, 26908, 50258,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 75698, 0, 3,
                                                                       74249, 49313, 74312,
                                                                       26908, 26968, 50348,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 75824, 0, 3,
                                                                       74312, 49358, 74375,
                                                                       26968, 27028, 50438,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 75950, 0, 3,
                                                                       74375, 49403, 74438,
                                                                       27028, 27088, 50528,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 76076, 0, 3,
                                                                       74438, 49448, 74501,
                                                                       27088, 27148, 50618,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 76202, 0, 3,
                                                                       74501, 49493, 74564,
                                                                       27148, 27208, 50708,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 76328, 0, 3,
                                                                       74564, 49538, 74627,
                                                                       27208, 27268, 50798,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 76454, 0, 3,
                                                                       74690, 49628, 74753,
                                                                       27388, 27448, 50888,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 76580, 0, 3,
                                                                       74753, 49673, 74816,
                                                                       27448, 27508, 50978,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 76706, 0, 3,
                                                                       74816, 49718, 74879,
                                                                       27508, 27568, 51068,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 76832, 0, 3,
                                                                       74879, 49763, 74942,
                                                                       27568, 27628, 51158,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 76958, 0, 3,
                                                                       74942, 49808, 75005,
                                                                       27628, 27688, 51248,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 77084, 0, 3,
                                                                       75005, 49853, 75068,
                                                                       27688, 27748, 51338,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 77210, 0, 3,
                                                                       75068, 49898, 75131,
                                                                       27748, 27808, 51428,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 77336, 0, 3,
                                                                       75131, 49943, 75194,
                                                                       27808, 27868, 51518,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 77462, 0, 3,
                                                                       75194, 49988, 75257,
                                                                       27868, 27928, 51608,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 77588, 0, 3,
                                                                       75320, 50078, 75446,
                                                                       28048, 28148, 51698,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 77798, 0, 3,
                                                                       75446, 50168, 75572,
                                                                       28148, 28248, 51848,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 78008, 0, 3,
                                                                       75572, 50258, 75698,
                                                                       28248, 28348, 51998,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 78218, 0, 3,
                                                                       75698, 50348, 75824,
                                                                       28348, 28448, 52148,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 78428, 0, 3,
                                                                       75824, 50438, 75950,
                                                                       28448, 28548, 52298,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 78638, 0, 3,
                                                                       75950, 50528, 76076,
                                                                       28548, 28648, 52448,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 78848, 0, 3,
                                                                       76076, 50618, 76202,
                                                                       28648, 28748, 52598,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 79058, 0, 3,
                                                                       76202, 50708, 76328,
                                                                       28748, 28848, 52748,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 79268, 0, 3,
                                                                       76454, 50888, 76580,
                                                                       29048, 29148, 52898,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 79478, 0, 3,
                                                                       76580, 50978, 76706,
                                                                       29148, 29248, 53048,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 79688, 0, 3,
                                                                       76706, 51068, 76832,
                                                                       29248, 29348, 53198,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 79898, 0, 3,
                                                                       76832, 51158, 76958,
                                                                       29348, 29448, 53348,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 80108, 0, 3,
                                                                       76958, 51248, 77084,
                                                                       29448, 29548, 53498,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 80318, 0, 3,
                                                                       77084, 51338, 77210,
                                                                       29548, 29648, 53648,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 80528, 0, 3,
                                                                       77210, 51428, 77336,
                                                                       29648, 29748, 53798,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 80738, 0, 3,
                                                                       77336, 51518, 77462,
                                                                       29748, 29848, 53948,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 80948, 0, 3,
                                                                       77588, 51698, 77798,
                                                                       30048, 30198, 54098,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 81263, 0, 3,
                                                                       77798, 51848, 78008,
                                                                       30198, 30348, 54323,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 81578, 0, 3,
                                                                       78008, 51998, 78218,
                                                                       30348, 30498, 54548,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 81893, 0, 3,
                                                                       78218, 52148, 78428,
                                                                       30498, 30648, 54773,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 82208, 0, 3,
                                                                       78428, 52298, 78638,
                                                                       30648, 30798, 54998,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 82523, 0, 3,
                                                                       78638, 52448, 78848,
                                                                       30798, 30948, 55223,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 82838, 0, 3,
                                                                       78848, 52598, 79058,
                                                                       30948, 31098, 55448,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 83153, 0, 3,
                                                                       79268, 52898, 79478,
                                                                       31398, 31548, 55673,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 83468, 0, 3,
                                                                       79478, 53048, 79688,
                                                                       31548, 31698, 55898,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 83783, 0, 3,
                                                                       79688, 53198, 79898,
                                                                       31698, 31848, 56123,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 84098, 0, 3,
                                                                       79898, 53348, 80108,
                                                                       31848, 31998, 56348,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 84413, 0, 3,
                                                                       80108, 53498, 80318,
                                                                       31998, 32148, 56573,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 84728, 0, 3,
                                                                       80318, 53648, 80528,
                                                                       32148, 32298, 56798,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 85043, 0, 3,
                                                                       80528, 53798, 80738,
                                                                       32298, 32448, 57023,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 85358, 0, 3,
                                                                       80948, 54098, 81263,
                                                                       32748, 32958, 57248,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 85799, 0, 3,
                                                                       81263, 54323, 81578,
                                                                       32958, 33168, 57563,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 86240, 0, 3,
                                                                       81578, 54548, 81893,
                                                                       33168, 33378, 57878,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 86681, 0, 3,
                                                                       81893, 54773, 82208,
                                                                       33378, 33588, 58193,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 87122, 0, 3,
                                                                       82208, 54998, 82523,
                                                                       33588, 33798, 58508,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 87563, 0, 3,
                                                                       82523, 55223, 82838,
                                                                       33798, 34008, 58823,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 88004, 0, 3,
                                                                       83153, 55673, 83468,
                                                                       34428, 34638, 59138,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 88445, 0, 3,
                                                                       83468, 55898, 83783,
                                                                       34638, 34848, 59453,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 88886, 0, 3,
                                                                       83783, 56123, 84098,
                                                                       34848, 35058, 59768,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 89327, 0, 3,
                                                                       84098, 56348, 84413,
                                                                       35058, 35268, 60083,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 89768, 0, 3,
                                                                       84413, 56573, 84728,
                                                                       35268, 35478, 60398,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 90209, 0, 3,
                                                                       84728, 56798, 85043,
                                                                       35478, 35688, 60713,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 90650, 0, 3,
                                                                       85358, 57248, 85799,
                                                                       36108, 36388, 61028,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 91238, 0, 3,
                                                                       85799, 57563, 86240,
                                                                       36388, 36668, 61448,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 91826, 0, 3,
                                                                       86240, 57878, 86681,
                                                                       36668, 36948, 61868,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 92414, 0, 3,
                                                                       86681, 58193, 87122,
                                                                       36948, 37228, 62288,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 93002, 0, 3,
                                                                       87122, 58508, 87563,
                                                                       37228, 37508, 62708,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 93590, 0, 3,
                                                                       88004, 59138, 88445,
                                                                       38068, 38348, 63128,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 94178, 0, 3,
                                                                       88445, 59453, 88886,
                                                                       38348, 38628, 63548,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 94766, 0, 3,
                                                                       88886, 59768, 89327,
                                                                       38628, 38908, 63968,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 95354, 0, 3,
                                                                       89327, 60083, 89768,
                                                                       38908, 39188, 64388,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 95942, 0, 3,
                                                                       89768, 60398, 90209,
                                                                       39188, 39468, 64808,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 96530, 0, 3,
                                                                       90650, 61028, 91238,
                                                                       40028, 40388, 65228,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 97286, 0, 3,
                                                                       91238, 61448, 91826,
                                                                       40388, 40748, 65768,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 98042, 0, 3,
                                                                       91826, 61868, 92414,
                                                                       40748, 41108, 66308,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 98798, 0, 3,
                                                                       92414, 62288, 93002,
                                                                       41108, 41468, 66848,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 99554, 0, 3,
                                                                       93590, 63128, 94178,
                                                                       42188, 42548, 67388,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 100310, 0, 3,
                                                                       94178, 63548, 94766,
                                                                       42548, 42908, 67928,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 101066, 0, 3,
                                                                       94766, 63968, 95354,
                                                                       42908, 43268, 68468,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 101822, 0, 3,
                                                                       95354, 64388, 95942,
                                                                       43268, 43628, 69008,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 102578, 0, 3,
                                                                       96530, 65228, 97286,
                                                                       44348, 44798, 69548,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 103523, 0, 3,
                                                                       97286, 65768, 98042,
                                                                       44798, 45248, 70223,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 104468, 0, 3,
                                                                       98042, 66308, 98798,
                                                                       45248, 45698, 70898,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 105413, 0, 3,
                                                                       99554, 67388, 100310,
                                                                       46598, 47048, 71573,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 106358, 0, 3,
                                                                       100310, 67928, 101066,
                                                                       47048, 47498, 72248,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 107303, 0, 3,
                                                                       101066, 68468, 101822,
                                                                       47498, 47948, 72923,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108248, 3, 48848,
                                                                       48863, 73640, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108276, 3, 48863,
                                                                       48878, 73661, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108304, 3, 48878,
                                                                       48893, 73682, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108332, 3, 48893,
                                                                       48908, 73703, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108360, 3, 48908,
                                                                       48923, 73724, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108388, 3, 48923,
                                                                       48938, 73745, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108416, 3, 48938,
                                                                       48953, 73766, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108444, 3, 48953,
                                                                       48968, 73787, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108472, 3, 48968,
                                                                       48983, 73808, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108500, 3, 49013,
                                                                       49028, 73871, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108528, 3, 49028,
                                                                       49043, 73892, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108556, 3, 49043,
                                                                       49058, 73913, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108584, 3, 49058,
                                                                       49073, 73934, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108612, 3, 49073,
                                                                       49088, 73955, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108640, 3, 49088,
                                                                       49103, 73976, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108668, 3, 49103,
                                                                       49118, 73997, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108696, 3, 49118,
                                                                       49133, 74018, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 108724, 3, 49133,
                                                                       49148, 74039, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 108752, 0, 3,
                                                                       108248, 73640, 108276,
                                                                       49178, 49223, 74186,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 108836, 0, 3,
                                                                       108276, 73661, 108304,
                                                                       49223, 49268, 74249,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 108920, 0, 3,
                                                                       108304, 73682, 108332,
                                                                       49268, 49313, 74312,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 109004, 0, 3,
                                                                       108332, 73703, 108360,
                                                                       49313, 49358, 74375,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 109088, 0, 3,
                                                                       108360, 73724, 108388,
                                                                       49358, 49403, 74438,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 109172, 0, 3,
                                                                       108388, 73745, 108416,
                                                                       49403, 49448, 74501,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 109256, 0, 3,
                                                                       108416, 73766, 108444,
                                                                       49448, 49493, 74564,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 109340, 0, 3,
                                                                       108444, 73787, 108472,
                                                                       49493, 49538, 74627,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 109424, 0, 3,
                                                                       108500, 73871, 108528,
                                                                       49628, 49673, 74816,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 109508, 0, 3,
                                                                       108528, 73892, 108556,
                                                                       49673, 49718, 74879,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 109592, 0, 3,
                                                                       108556, 73913, 108584,
                                                                       49718, 49763, 74942,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 109676, 0, 3,
                                                                       108584, 73934, 108612,
                                                                       49763, 49808, 75005,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 109760, 0, 3,
                                                                       108612, 73955, 108640,
                                                                       49808, 49853, 75068,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 109844, 0, 3,
                                                                       108640, 73976, 108668,
                                                                       49853, 49898, 75131,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 109928, 0, 3,
                                                                       108668, 73997, 108696,
                                                                       49898, 49943, 75194,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 110012, 0, 3,
                                                                       108696, 74018, 108724,
                                                                       49943, 49988, 75257,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 110096, 0, 3,
                                                                       108752, 74186, 108836,
                                                                       50078, 50168, 75572,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 110264, 0, 3,
                                                                       108836, 74249, 108920,
                                                                       50168, 50258, 75698,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 110432, 0, 3,
                                                                       108920, 74312, 109004,
                                                                       50258, 50348, 75824,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 110600, 0, 3,
                                                                       109004, 74375, 109088,
                                                                       50348, 50438, 75950,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 110768, 0, 3,
                                                                       109088, 74438, 109172,
                                                                       50438, 50528, 76076,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 110936, 0, 3,
                                                                       109172, 74501, 109256,
                                                                       50528, 50618, 76202,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 111104, 0, 3,
                                                                       109256, 74564, 109340,
                                                                       50618, 50708, 76328,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 111272, 0, 3,
                                                                       109424, 74816, 109508,
                                                                       50888, 50978, 76706,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 111440, 0, 3,
                                                                       109508, 74879, 109592,
                                                                       50978, 51068, 76832,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 111608, 0, 3,
                                                                       109592, 74942, 109676,
                                                                       51068, 51158, 76958,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 111776, 0, 3,
                                                                       109676, 75005, 109760,
                                                                       51158, 51248, 77084,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 111944, 0, 3,
                                                                       109760, 75068, 109844,
                                                                       51248, 51338, 77210,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 112112, 0, 3,
                                                                       109844, 75131, 109928,
                                                                       51338, 51428, 77336,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 112280, 0, 3,
                                                                       109928, 75194, 110012,
                                                                       51428, 51518, 77462,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 112448, 0, 3,
                                                                       110096, 75572, 110264,
                                                                       51698, 51848, 78008,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 112728, 0, 3,
                                                                       110264, 75698, 110432,
                                                                       51848, 51998, 78218,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 113008, 0, 3,
                                                                       110432, 75824, 110600,
                                                                       51998, 52148, 78428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 113288, 0, 3,
                                                                       110600, 75950, 110768,
                                                                       52148, 52298, 78638,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 113568, 0, 3,
                                                                       110768, 76076, 110936,
                                                                       52298, 52448, 78848,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 113848, 0, 3,
                                                                       110936, 76202, 111104,
                                                                       52448, 52598, 79058,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 114128, 0, 3,
                                                                       111272, 76706, 111440,
                                                                       52898, 53048, 79688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 114408, 0, 3,
                                                                       111440, 76832, 111608,
                                                                       53048, 53198, 79898,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 114688, 0, 3,
                                                                       111608, 76958, 111776,
                                                                       53198, 53348, 80108,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 114968, 0, 3,
                                                                       111776, 77084, 111944,
                                                                       53348, 53498, 80318,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 115248, 0, 3,
                                                                       111944, 77210, 112112,
                                                                       53498, 53648, 80528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 115528, 0, 3,
                                                                       112112, 77336, 112280,
                                                                       53648, 53798, 80738,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 115808, 0, 3,
                                                                       112448, 78008, 112728,
                                                                       54098, 54323, 81578,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 116228, 0, 3,
                                                                       112728, 78218, 113008,
                                                                       54323, 54548, 81893,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 116648, 0, 3,
                                                                       113008, 78428, 113288,
                                                                       54548, 54773, 82208,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 117068, 0, 3,
                                                                       113288, 78638, 113568,
                                                                       54773, 54998, 82523,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 117488, 0, 3,
                                                                       113568, 78848, 113848,
                                                                       54998, 55223, 82838,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 117908, 0, 3,
                                                                       114128, 79688, 114408,
                                                                       55673, 55898, 83783,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 118328, 0, 3,
                                                                       114408, 79898, 114688,
                                                                       55898, 56123, 84098,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 118748, 0, 3,
                                                                       114688, 80108, 114968,
                                                                       56123, 56348, 84413,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 119168, 0, 3,
                                                                       114968, 80318, 115248,
                                                                       56348, 56573, 84728,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 119588, 0, 3,
                                                                       115248, 80528, 115528,
                                                                       56573, 56798, 85043,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 120008, 0, 3,
                                                                       115808, 81578, 116228,
                                                                       57248, 57563, 86240,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 120596, 0, 3,
                                                                       116228, 81893, 116648,
                                                                       57563, 57878, 86681,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 121184, 0, 3,
                                                                       116648, 82208, 117068,
                                                                       57878, 58193, 87122,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 121772, 0, 3,
                                                                       117068, 82523, 117488,
                                                                       58193, 58508, 87563,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 122360, 0, 3,
                                                                       117908, 83783, 118328,
                                                                       59138, 59453, 88886,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 122948, 0, 3,
                                                                       118328, 84098, 118748,
                                                                       59453, 59768, 89327,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 123536, 0, 3,
                                                                       118748, 84413, 119168,
                                                                       59768, 60083, 89768,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 124124, 0, 3,
                                                                       119168, 84728, 119588,
                                                                       60083, 60398, 90209,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 124712, 0, 3,
                                                                       120008, 86240, 120596,
                                                                       61028, 61448, 91826,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 125496, 0, 3,
                                                                       120596, 86681, 121184,
                                                                       61448, 61868, 92414,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 126280, 0, 3,
                                                                       121184, 87122, 121772,
                                                                       61868, 62288, 93002,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 127064, 0, 3,
                                                                       122360, 88886, 122948,
                                                                       63128, 63548, 94766,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 127848, 0, 3,
                                                                       122948, 89327, 123536,
                                                                       63548, 63968, 95354,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 128632, 0, 3,
                                                                       123536, 89768, 124124,
                                                                       63968, 64388, 95942,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 129416, 0, 3,
                                                                       124712, 91826, 125496,
                                                                       65228, 65768, 98042,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 130424, 0, 3,
                                                                       125496, 92414, 126280,
                                                                       65768, 66308, 98798,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 131432, 0, 3,
                                                                       127064, 94766, 127848,
                                                                       67388, 67928, 101066,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 132440, 0, 3,
                                                                       127848, 95354, 128632,
                                                                       67928, 68468, 101822,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 133448, 0, 3,
                                                                       129416, 98042, 130424,
                                                                       69548, 70223, 104468,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 134708, 0, 3,
                                                                       131432, 101066, 132440,
                                                                       71573, 72248, 107303,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 135968, 3, 73598,
                                                                       73619, 108248, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136004, 3, 73619,
                                                                       73640, 108276, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136040, 3, 73640,
                                                                       73661, 108304, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136076, 3, 73661,
                                                                       73682, 108332, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136112, 3, 73682,
                                                                       73703, 108360, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136148, 3, 73703,
                                                                       73724, 108388, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136184, 3, 73724,
                                                                       73745, 108416, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136220, 3, 73745,
                                                                       73766, 108444, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136256, 3, 73766,
                                                                       73787, 108472, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136292, 3, 73829,
                                                                       73850, 108500, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136328, 3, 73850,
                                                                       73871, 108528, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136364, 3, 73871,
                                                                       73892, 108556, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136400, 3, 73892,
                                                                       73913, 108584, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136436, 3, 73913,
                                                                       73934, 108612, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136472, 3, 73934,
                                                                       73955, 108640, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136508, 3, 73955,
                                                                       73976, 108668, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136544, 3, 73976,
                                                                       73997, 108696, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 136580, 3, 73997,
                                                                       74018, 108724, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 136616, 0, 3,
                                                                       135968, 108248, 136004,
                                                                       74060, 74123, 108752,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 136724, 0, 3,
                                                                       136004, 108276, 136040,
                                                                       74123, 74186, 108836,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 136832, 0, 3,
                                                                       136040, 108304, 136076,
                                                                       74186, 74249, 108920,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 136940, 0, 3,
                                                                       136076, 108332, 136112,
                                                                       74249, 74312, 109004,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 137048, 0, 3,
                                                                       136112, 108360, 136148,
                                                                       74312, 74375, 109088,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 137156, 0, 3,
                                                                       136148, 108388, 136184,
                                                                       74375, 74438, 109172,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 137264, 0, 3,
                                                                       136184, 108416, 136220,
                                                                       74438, 74501, 109256,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 137372, 0, 3,
                                                                       136220, 108444, 136256,
                                                                       74501, 74564, 109340,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 137480, 0, 3,
                                                                       136292, 108500, 136328,
                                                                       74690, 74753, 109424,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 137588, 0, 3,
                                                                       136328, 108528, 136364,
                                                                       74753, 74816, 109508,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 137696, 0, 3,
                                                                       136364, 108556, 136400,
                                                                       74816, 74879, 109592,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 137804, 0, 3,
                                                                       136400, 108584, 136436,
                                                                       74879, 74942, 109676,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 137912, 0, 3,
                                                                       136436, 108612, 136472,
                                                                       74942, 75005, 109760,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 138020, 0, 3,
                                                                       136472, 108640, 136508,
                                                                       75005, 75068, 109844,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 138128, 0, 3,
                                                                       136508, 108668, 136544,
                                                                       75068, 75131, 109928,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 138236, 0, 3,
                                                                       136544, 108696, 136580,
                                                                       75131, 75194, 110012,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 138344, 0, 3,
                                                                       136616, 108752, 136724,
                                                                       75320, 75446, 110096,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 138560, 0, 3,
                                                                       136724, 108836, 136832,
                                                                       75446, 75572, 110264,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 138776, 0, 3,
                                                                       136832, 108920, 136940,
                                                                       75572, 75698, 110432,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 138992, 0, 3,
                                                                       136940, 109004, 137048,
                                                                       75698, 75824, 110600,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 139208, 0, 3,
                                                                       137048, 109088, 137156,
                                                                       75824, 75950, 110768,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 139424, 0, 3,
                                                                       137156, 109172, 137264,
                                                                       75950, 76076, 110936,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 139640, 0, 3,
                                                                       137264, 109256, 137372,
                                                                       76076, 76202, 111104,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 139856, 0, 3,
                                                                       137480, 109424, 137588,
                                                                       76454, 76580, 111272,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 140072, 0, 3,
                                                                       137588, 109508, 137696,
                                                                       76580, 76706, 111440,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 140288, 0, 3,
                                                                       137696, 109592, 137804,
                                                                       76706, 76832, 111608,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 140504, 0, 3,
                                                                       137804, 109676, 137912,
                                                                       76832, 76958, 111776,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 140720, 0, 3,
                                                                       137912, 109760, 138020,
                                                                       76958, 77084, 111944,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 140936, 0, 3,
                                                                       138020, 109844, 138128,
                                                                       77084, 77210, 112112,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 141152, 0, 3,
                                                                       138128, 109928, 138236,
                                                                       77210, 77336, 112280,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 141368, 0, 3,
                                                                       138344, 110096, 138560,
                                                                       77588, 77798, 112448,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 141728, 0, 3,
                                                                       138560, 110264, 138776,
                                                                       77798, 78008, 112728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 142088, 0, 3,
                                                                       138776, 110432, 138992,
                                                                       78008, 78218, 113008,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 142448, 0, 3,
                                                                       138992, 110600, 139208,
                                                                       78218, 78428, 113288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 142808, 0, 3,
                                                                       139208, 110768, 139424,
                                                                       78428, 78638, 113568,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 143168, 0, 3,
                                                                       139424, 110936, 139640,
                                                                       78638, 78848, 113848,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 143528, 0, 3,
                                                                       139856, 111272, 140072,
                                                                       79268, 79478, 114128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 143888, 0, 3,
                                                                       140072, 111440, 140288,
                                                                       79478, 79688, 114408,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 144248, 0, 3,
                                                                       140288, 111608, 140504,
                                                                       79688, 79898, 114688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 144608, 0, 3,
                                                                       140504, 111776, 140720,
                                                                       79898, 80108, 114968,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 144968, 0, 3,
                                                                       140720, 111944, 140936,
                                                                       80108, 80318, 115248,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 145328, 0, 3,
                                                                       140936, 112112, 141152,
                                                                       80318, 80528, 115528,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 145688, 0, 3,
                                                                       141368, 112448, 141728,
                                                                       80948, 81263, 115808,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 146228, 0, 3,
                                                                       141728, 112728, 142088,
                                                                       81263, 81578, 116228,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 146768, 0, 3,
                                                                       142088, 113008, 142448,
                                                                       81578, 81893, 116648,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 147308, 0, 3,
                                                                       142448, 113288, 142808,
                                                                       81893, 82208, 117068,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 147848, 0, 3,
                                                                       142808, 113568, 143168,
                                                                       82208, 82523, 117488,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 148388, 0, 3,
                                                                       143528, 114128, 143888,
                                                                       83153, 83468, 117908,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 148928, 0, 3,
                                                                       143888, 114408, 144248,
                                                                       83468, 83783, 118328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 149468, 0, 3,
                                                                       144248, 114688, 144608,
                                                                       83783, 84098, 118748,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 150008, 0, 3,
                                                                       144608, 114968, 144968,
                                                                       84098, 84413, 119168,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 150548, 0, 3,
                                                                       144968, 115248, 145328,
                                                                       84413, 84728, 119588,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 151088, 0, 3,
                                                                       145688, 115808, 146228,
                                                                       85358, 85799, 120008,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 151844, 0, 3,
                                                                       146228, 116228, 146768,
                                                                       85799, 86240, 120596,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 152600, 0, 3,
                                                                       146768, 116648, 147308,
                                                                       86240, 86681, 121184,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 153356, 0, 3,
                                                                       147308, 117068, 147848,
                                                                       86681, 87122, 121772,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 154112, 0, 3,
                                                                       148388, 117908, 148928,
                                                                       88004, 88445, 122360,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 154868, 0, 3,
                                                                       148928, 118328, 149468,
                                                                       88445, 88886, 122948,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 155624, 0, 3,
                                                                       149468, 118748, 150008,
                                                                       88886, 89327, 123536,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 156380, 0, 3,
                                                                       150008, 119168, 150548,
                                                                       89327, 89768, 124124,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 157136, 0, 3,
                                                                       151088, 120008, 151844,
                                                                       90650, 91238, 124712,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 158144, 0, 3,
                                                                       151844, 120596, 152600,
                                                                       91238, 91826, 125496,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 159152, 0, 3,
                                                                       152600, 121184, 153356,
                                                                       91826, 92414, 126280,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 160160, 0, 3,
                                                                       154112, 122360, 154868,
                                                                       93590, 94178, 127064,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 161168, 0, 3,
                                                                       154868, 122948, 155624,
                                                                       94178, 94766, 127848,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 162176, 0, 3,
                                                                       155624, 123536, 156380,
                                                                       94766, 95354, 128632,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 163184, 0, 3,
                                                                       157136, 124712, 158144,
                                                                       96530, 97286, 129416,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 164480, 0, 3,
                                                                       158144, 125496, 159152,
                                                                       97286, 98042, 130424,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 165776, 0, 3,
                                                                       160160, 127064, 161168,
                                                                       99554, 100310, 131432,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 167072, 0, 3,
                                                                       161168, 127848, 162176,
                                                                       100310, 101066, 132440,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 168368, 0, 3,
                                                                       163184, 129416, 164480,
                                                                       102578, 103523, 133448,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 169988, 0, 3,
                                                                       165776, 131432, 167072,
                                                                       105413, 106358, 134708,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 171608, 151088, 756, ncols);

                    simdfunc::contract_primitives(buffer, 172679, 154112, 756, ncols);

                    simdfunc::contract_primitives(buffer, 173750, 157136, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 175178, 160160, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 176606, 163184, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 178442, 165776, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 180278, 168368, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 182573, 169988, 1620, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 172364, 171608, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 173435, 172679, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 174758, 173750, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 176186, 175178, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 177902, 176606, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 179738, 178442, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 181898, 180278, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 184193, 182573, 45, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 184868, 172364, 174758, 15,
                                             nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 185813, 173435, 176186, 15,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 186758, 174758, 177902, 15,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 188018, 176186, 179738, 15,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 189278, 177902, 181898, 15,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 190898, 179738, 184193, 15,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 192518, 184868, 186758, 15,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 194408, 185813, 188018, 15,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 196298, 186758, 189278, 15,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 198818, 188018, 190898, 15,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 201338, 192518, 196298, 15,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 204488, 194408, 198818, 15,
                                             nmax);

        simdtrf::transform_f_inner(buffer, 207638, 204488, 21, 15, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 207638, 105, nmax);

        simdtrf::transform_f_inner(buffer, 207638, 201338, 21, 15, nmax);

        simdtrf::transform_h_outer(values + 1155 * nvalues + n * npairs, nvalues, buffer, 207638,
                                   105, nmax);
    }

    for (size_t m = 0; m < 2310; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
