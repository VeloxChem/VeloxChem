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


#include "SimdThreeCenterElectronRepulsionRsRecIGG.hpp"

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
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKF.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLD.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransferMP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_igg_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_igg_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 135678, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2106 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 135678, 82662, 10446, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 14,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 22, 3, 14,
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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4298, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4301, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4304, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4307, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4310, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4313, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4316, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4319, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4322, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4325, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4328, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4331, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4334, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4337, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4340, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4343, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4346, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4349, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4352, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4355, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4358, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4361, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4364, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4367, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4370, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4373, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4376, 3, 9, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4385, 3, 10, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4394, 3, 11, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4403, 3, 12, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4412, 3, 13, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4421, 3, 14, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4430, 3, 15, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4439, 3, 16, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4448, 3, 17, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4457, 3, 18, 71,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4466, 3, 19, 74,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4475, 3, 20, 77,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4484, 3, 25, 86,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4493, 3, 26, 89,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4502, 3, 27, 92,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4511, 3, 28, 95,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4520, 3, 29, 98,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4529, 3, 30, 101,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4538, 3, 31, 104,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4547, 3, 32, 107,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4556, 3, 33, 110,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4565, 3, 34, 113,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4574, 3, 35, 116,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4583, 3, 36, 119,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4592, 3, 44, 134,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4610, 3, 47, 140,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4628, 3, 50, 146,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4646, 3, 53, 152,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4664, 3, 56, 158,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4682, 3, 59, 164,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4700, 3, 62, 170,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4718, 3, 65, 176,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4736, 3, 68, 182,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4754, 3, 71, 188,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4772, 3, 74, 194,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4790, 3, 86, 212,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4808, 3, 89, 218,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4826, 3, 92, 224,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4844, 3, 95, 230,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4862, 3, 98, 236,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4880, 3, 101, 242,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4898, 3, 104, 248,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4916, 3, 107, 254,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4934, 3, 110, 260,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4952, 3, 113, 266,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4970, 3, 116, 272,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4988, 3, 134, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5018, 3, 140, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5048, 3, 146, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5078, 3, 152, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5108, 3, 158, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5138, 3, 164, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5168, 3, 170, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5198, 3, 176, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5228, 3, 182, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5258, 3, 188, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5288, 3, 212, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5318, 3, 218, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5348, 3, 224, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5378, 3, 230, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5408, 3, 236, 458,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5438, 3, 242, 468,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5468, 3, 248, 478,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5498, 3, 254, 488,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5528, 3, 260, 498,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5558, 3, 266, 508,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5588, 3, 298, 548,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5633, 3, 308, 563,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5678, 3, 318, 578,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5723, 3, 328, 593,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5768, 3, 338, 608,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5813, 3, 348, 623,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5858, 3, 358, 638,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5903, 3, 368, 653,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5948, 3, 378, 668,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5993, 3, 418, 713,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6038, 3, 428, 728,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6083, 3, 438, 743,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6128, 3, 448, 758,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6173, 3, 458, 773,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6218, 3, 468, 788,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6263, 3, 478, 803,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6308, 3, 488, 818,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6353, 3, 498, 833,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6398, 3, 548, 890,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6461, 3, 563, 911,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6524, 3, 578, 932,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6587, 3, 593, 953,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6650, 3, 608, 974,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6713, 3, 623, 995,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6776, 3, 638,
                                                                       1016, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6839, 3, 653,
                                                                       1037, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6902, 3, 713,
                                                                       1100, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6965, 3, 728,
                                                                       1121, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7028, 3, 743,
                                                                       1142, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7091, 3, 758,
                                                                       1163, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7154, 3, 773,
                                                                       1184, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7217, 3, 788,
                                                                       1205, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7280, 3, 803,
                                                                       1226, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7343, 3, 818,
                                                                       1247, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7406, 3, 890,
                                                                       1324, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7490, 3, 911,
                                                                       1352, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7574, 3, 932,
                                                                       1380, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7658, 3, 953,
                                                                       1408, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7742, 3, 974,
                                                                       1436, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7826, 3, 995,
                                                                       1464, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7910, 3, 1016,
                                                                       1492, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7994, 3, 1100,
                                                                       1576, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8078, 3, 1121,
                                                                       1604, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8162, 3, 1142,
                                                                       1632, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8246, 3, 1163,
                                                                       1660, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8330, 3, 1184,
                                                                       1688, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8414, 3, 1205,
                                                                       1716, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8498, 3, 1226,
                                                                       1744, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8582, 3, 1324,
                                                                       1844, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8690, 3, 1352,
                                                                       1880, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8798, 3, 1380,
                                                                       1916, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8906, 3, 1408,
                                                                       1952, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9014, 3, 1436,
                                                                       1988, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9122, 3, 1464,
                                                                       2024, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9230, 3, 1576,
                                                                       2132, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9338, 3, 1604,
                                                                       2168, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9446, 3, 1632,
                                                                       2204, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9554, 3, 1660,
                                                                       2240, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9662, 3, 1688,
                                                                       2276, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9770, 3, 1716,
                                                                       2312, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9878, 3, 1844,
                                                                       2438, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10013, 3, 1880,
                                                                       2483, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10148, 3, 1916,
                                                                       2528, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10283, 3, 1952,
                                                                       2573, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10418, 3, 1988,
                                                                       2618, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10553, 3, 2132,
                                                                       2753, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10688, 3, 2168,
                                                                       2798, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10823, 3, 2204,
                                                                       2843, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10958, 3, 2240,
                                                                       2888, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11093, 3, 2276,
                                                                       2933, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11228, 3, 2438,
                                                                       3088, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11393, 3, 2483,
                                                                       3143, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11558, 3, 2528,
                                                                       3198, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11723, 3, 2573,
                                                                       3253, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11888, 3, 2753,
                                                                       3418, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12053, 3, 2798,
                                                                       3473, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12218, 3, 2843,
                                                                       3528, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12383, 3, 2888,
                                                                       3583, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 12548, 3, 3088,
                                                                       3770, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 12746, 3, 3143,
                                                                       3836, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 12944, 3, 3198,
                                                                       3902, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13142, 3, 3418,
                                                                       4100, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13340, 3, 3473,
                                                                       4166, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13538, 3, 3528,
                                                                       4232, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13736, 3, 7, 8,
                                                                       4298, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13742, 3, 8, 9,
                                                                       4301, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13748, 3, 9, 10,
                                                                       4304, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13754, 3, 10, 11,
                                                                       4307, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13760, 3, 11, 12,
                                                                       4310, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13766, 3, 12, 13,
                                                                       4313, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13772, 3, 13, 14,
                                                                       4316, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13778, 3, 14, 15,
                                                                       4319, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13784, 3, 15, 16,
                                                                       4322, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13790, 3, 16, 17,
                                                                       4325, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13796, 3, 17, 18,
                                                                       4328, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13802, 3, 18, 19,
                                                                       4331, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13808, 3, 19, 20,
                                                                       4334, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13814, 3, 23, 24,
                                                                       4337, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13820, 3, 24, 25,
                                                                       4340, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13826, 3, 25, 26,
                                                                       4343, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13832, 3, 26, 27,
                                                                       4346, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13838, 3, 27, 28,
                                                                       4349, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13844, 3, 28, 29,
                                                                       4352, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13850, 3, 29, 30,
                                                                       4355, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13856, 3, 30, 31,
                                                                       4358, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13862, 3, 31, 32,
                                                                       4361, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13868, 3, 32, 33,
                                                                       4364, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13874, 3, 33, 34,
                                                                       4367, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13880, 3, 34, 35,
                                                                       4370, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13886, 3, 35, 36,
                                                                       4373, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13892, 0, 3,
                                                                       13736, 4298, 13742, 4376,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13910, 0, 3,
                                                                       13742, 4301, 13748, 4385,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13928, 0, 3,
                                                                       13748, 4304, 13754, 4394,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13946, 0, 3,
                                                                       13754, 4307, 13760, 4403,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13964, 0, 3,
                                                                       13760, 4310, 13766, 4412,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 13982, 0, 3,
                                                                       13766, 4313, 13772, 4421,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14000, 0, 3,
                                                                       13772, 4316, 13778, 4430,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14018, 0, 3,
                                                                       13778, 4319, 13784, 4439,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14036, 0, 3,
                                                                       13784, 4322, 13790, 4448,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14054, 0, 3,
                                                                       13790, 4325, 13796, 4457,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14072, 0, 3,
                                                                       13796, 4328, 13802, 4466,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14090, 0, 3,
                                                                       13802, 4331, 13808, 4475,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14108, 0, 3,
                                                                       13814, 4337, 13820, 4484,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14126, 0, 3,
                                                                       13820, 4340, 13826, 4493,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14144, 0, 3,
                                                                       13826, 4343, 13832, 4502,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14162, 0, 3,
                                                                       13832, 4346, 13838, 4511,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14180, 0, 3,
                                                                       13838, 4349, 13844, 4520,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14198, 0, 3,
                                                                       13844, 4352, 13850, 4529,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14216, 0, 3,
                                                                       13850, 4355, 13856, 4538,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14234, 0, 3,
                                                                       13856, 4358, 13862, 4547,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14252, 0, 3,
                                                                       13862, 4361, 13868, 4556,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14270, 0, 3,
                                                                       13868, 4364, 13874, 4565,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14288, 0, 3,
                                                                       13874, 4367, 13880, 4574,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 14306, 0, 3,
                                                                       13880, 4370, 13886, 4583,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14324, 0, 3,
                                                                       13892, 4376, 13910, 122,
                                                                       128, 4592, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14360, 0, 3,
                                                                       13910, 4385, 13928, 128,
                                                                       134, 4610, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14396, 0, 3,
                                                                       13928, 4394, 13946, 134,
                                                                       140, 4628, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14432, 0, 3,
                                                                       13946, 4403, 13964, 140,
                                                                       146, 4646, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14468, 0, 3,
                                                                       13964, 4412, 13982, 146,
                                                                       152, 4664, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14504, 0, 3,
                                                                       13982, 4421, 14000, 152,
                                                                       158, 4682, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14540, 0, 3,
                                                                       14000, 4430, 14018, 158,
                                                                       164, 4700, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14576, 0, 3,
                                                                       14018, 4439, 14036, 164,
                                                                       170, 4718, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14612, 0, 3,
                                                                       14036, 4448, 14054, 170,
                                                                       176, 4736, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14648, 0, 3,
                                                                       14054, 4457, 14072, 176,
                                                                       182, 4754, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14684, 0, 3,
                                                                       14072, 4466, 14090, 182,
                                                                       188, 4772, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14720, 0, 3,
                                                                       14108, 4484, 14126, 200,
                                                                       206, 4790, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14756, 0, 3,
                                                                       14126, 4493, 14144, 206,
                                                                       212, 4808, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14792, 0, 3,
                                                                       14144, 4502, 14162, 212,
                                                                       218, 4826, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14828, 0, 3,
                                                                       14162, 4511, 14180, 218,
                                                                       224, 4844, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14864, 0, 3,
                                                                       14180, 4520, 14198, 224,
                                                                       230, 4862, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14900, 0, 3,
                                                                       14198, 4529, 14216, 230,
                                                                       236, 4880, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14936, 0, 3,
                                                                       14216, 4538, 14234, 236,
                                                                       242, 4898, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 14972, 0, 3,
                                                                       14234, 4547, 14252, 242,
                                                                       248, 4916, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15008, 0, 3,
                                                                       14252, 4556, 14270, 248,
                                                                       254, 4934, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15044, 0, 3,
                                                                       14270, 4565, 14288, 254,
                                                                       260, 4952, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15080, 0, 3,
                                                                       14288, 4574, 14306, 260,
                                                                       266, 4970, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15116, 0, 3,
                                                                       14324, 4592, 14360, 278,
                                                                       288, 4988, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15176, 0, 3,
                                                                       14360, 4610, 14396, 288,
                                                                       298, 5018, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15236, 0, 3,
                                                                       14396, 4628, 14432, 298,
                                                                       308, 5048, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15296, 0, 3,
                                                                       14432, 4646, 14468, 308,
                                                                       318, 5078, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15356, 0, 3,
                                                                       14468, 4664, 14504, 318,
                                                                       328, 5108, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15416, 0, 3,
                                                                       14504, 4682, 14540, 328,
                                                                       338, 5138, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15476, 0, 3,
                                                                       14540, 4700, 14576, 338,
                                                                       348, 5168, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15536, 0, 3,
                                                                       14576, 4718, 14612, 348,
                                                                       358, 5198, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15596, 0, 3,
                                                                       14612, 4736, 14648, 358,
                                                                       368, 5228, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15656, 0, 3,
                                                                       14648, 4754, 14684, 368,
                                                                       378, 5258, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15716, 0, 3,
                                                                       14720, 4790, 14756, 398,
                                                                       408, 5288, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15776, 0, 3,
                                                                       14756, 4808, 14792, 408,
                                                                       418, 5318, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15836, 0, 3,
                                                                       14792, 4826, 14828, 418,
                                                                       428, 5348, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15896, 0, 3,
                                                                       14828, 4844, 14864, 428,
                                                                       438, 5378, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 15956, 0, 3,
                                                                       14864, 4862, 14900, 438,
                                                                       448, 5408, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16016, 0, 3,
                                                                       14900, 4880, 14936, 448,
                                                                       458, 5438, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16076, 0, 3,
                                                                       14936, 4898, 14972, 458,
                                                                       468, 5468, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16136, 0, 3,
                                                                       14972, 4916, 15008, 468,
                                                                       478, 5498, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16196, 0, 3,
                                                                       15008, 4934, 15044, 478,
                                                                       488, 5528, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16256, 0, 3,
                                                                       15044, 4952, 15080, 488,
                                                                       498, 5558, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16316, 0, 3,
                                                                       15116, 4988, 15176, 518,
                                                                       533, 5588, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16406, 0, 3,
                                                                       15176, 5018, 15236, 533,
                                                                       548, 5633, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16496, 0, 3,
                                                                       15236, 5048, 15296, 548,
                                                                       563, 5678, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16586, 0, 3,
                                                                       15296, 5078, 15356, 563,
                                                                       578, 5723, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16676, 0, 3,
                                                                       15356, 5108, 15416, 578,
                                                                       593, 5768, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16766, 0, 3,
                                                                       15416, 5138, 15476, 593,
                                                                       608, 5813, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16856, 0, 3,
                                                                       15476, 5168, 15536, 608,
                                                                       623, 5858, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 16946, 0, 3,
                                                                       15536, 5198, 15596, 623,
                                                                       638, 5903, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17036, 0, 3,
                                                                       15596, 5228, 15656, 638,
                                                                       653, 5948, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17126, 0, 3,
                                                                       15716, 5288, 15776, 683,
                                                                       698, 5993, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17216, 0, 3,
                                                                       15776, 5318, 15836, 698,
                                                                       713, 6038, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17306, 0, 3,
                                                                       15836, 5348, 15896, 713,
                                                                       728, 6083, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17396, 0, 3,
                                                                       15896, 5378, 15956, 728,
                                                                       743, 6128, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17486, 0, 3,
                                                                       15956, 5408, 16016, 743,
                                                                       758, 6173, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17576, 0, 3,
                                                                       16016, 5438, 16076, 758,
                                                                       773, 6218, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17666, 0, 3,
                                                                       16076, 5468, 16136, 773,
                                                                       788, 6263, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17756, 0, 3,
                                                                       16136, 5498, 16196, 788,
                                                                       803, 6308, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17846, 0, 3,
                                                                       16196, 5528, 16256, 803,
                                                                       818, 6353, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 17936, 0, 3,
                                                                       16316, 5588, 16406, 848,
                                                                       869, 6398, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18062, 0, 3,
                                                                       16406, 5633, 16496, 869,
                                                                       890, 6461, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18188, 0, 3,
                                                                       16496, 5678, 16586, 890,
                                                                       911, 6524, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18314, 0, 3,
                                                                       16586, 5723, 16676, 911,
                                                                       932, 6587, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18440, 0, 3,
                                                                       16676, 5768, 16766, 932,
                                                                       953, 6650, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18566, 0, 3,
                                                                       16766, 5813, 16856, 953,
                                                                       974, 6713, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18692, 0, 3,
                                                                       16856, 5858, 16946, 974,
                                                                       995, 6776, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18818, 0, 3,
                                                                       16946, 5903, 17036, 995,
                                                                       1016, 6839, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18944, 0, 3,
                                                                       17126, 5993, 17216, 1058,
                                                                       1079, 6902, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19070, 0, 3,
                                                                       17216, 6038, 17306, 1079,
                                                                       1100, 6965, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19196, 0, 3,
                                                                       17306, 6083, 17396, 1100,
                                                                       1121, 7028, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19322, 0, 3,
                                                                       17396, 6128, 17486, 1121,
                                                                       1142, 7091, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19448, 0, 3,
                                                                       17486, 6173, 17576, 1142,
                                                                       1163, 7154, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19574, 0, 3,
                                                                       17576, 6218, 17666, 1163,
                                                                       1184, 7217, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19700, 0, 3,
                                                                       17666, 6263, 17756, 1184,
                                                                       1205, 7280, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19826, 0, 3,
                                                                       17756, 6308, 17846, 1205,
                                                                       1226, 7343, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19952, 0, 3,
                                                                       17936, 6398, 18062, 1268,
                                                                       1296, 7406, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20120, 0, 3,
                                                                       18062, 6461, 18188, 1296,
                                                                       1324, 7490, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20288, 0, 3,
                                                                       18188, 6524, 18314, 1324,
                                                                       1352, 7574, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20456, 0, 3,
                                                                       18314, 6587, 18440, 1352,
                                                                       1380, 7658, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20624, 0, 3,
                                                                       18440, 6650, 18566, 1380,
                                                                       1408, 7742, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20792, 0, 3,
                                                                       18566, 6713, 18692, 1408,
                                                                       1436, 7826, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20960, 0, 3,
                                                                       18692, 6776, 18818, 1436,
                                                                       1464, 7910, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 21128, 0, 3,
                                                                       18944, 6902, 19070, 1520,
                                                                       1548, 7994, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 21296, 0, 3,
                                                                       19070, 6965, 19196, 1548,
                                                                       1576, 8078, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 21464, 0, 3,
                                                                       19196, 7028, 19322, 1576,
                                                                       1604, 8162, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 21632, 0, 3,
                                                                       19322, 7091, 19448, 1604,
                                                                       1632, 8246, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 21800, 0, 3,
                                                                       19448, 7154, 19574, 1632,
                                                                       1660, 8330, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 21968, 0, 3,
                                                                       19574, 7217, 19700, 1660,
                                                                       1688, 8414, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22136, 0, 3,
                                                                       19700, 7280, 19826, 1688,
                                                                       1716, 8498, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 22304, 0, 3,
                                                                       19952, 7406, 20120, 1772,
                                                                       1808, 8582, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 22520, 0, 3,
                                                                       20120, 7490, 20288, 1808,
                                                                       1844, 8690, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 22736, 0, 3,
                                                                       20288, 7574, 20456, 1844,
                                                                       1880, 8798, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 22952, 0, 3,
                                                                       20456, 7658, 20624, 1880,
                                                                       1916, 8906, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 23168, 0, 3,
                                                                       20624, 7742, 20792, 1916,
                                                                       1952, 9014, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 23384, 0, 3,
                                                                       20792, 7826, 20960, 1952,
                                                                       1988, 9122, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 23600, 0, 3,
                                                                       21128, 7994, 21296, 2060,
                                                                       2096, 9230, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 23816, 0, 3,
                                                                       21296, 8078, 21464, 2096,
                                                                       2132, 9338, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 24032, 0, 3,
                                                                       21464, 8162, 21632, 2132,
                                                                       2168, 9446, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 24248, 0, 3,
                                                                       21632, 8246, 21800, 2168,
                                                                       2204, 9554, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 24464, 0, 3,
                                                                       21800, 8330, 21968, 2204,
                                                                       2240, 9662, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 24680, 0, 3,
                                                                       21968, 8414, 22136, 2240,
                                                                       2276, 9770, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 24896, 0, 3,
                                                                       22304, 8582, 22520, 2348,
                                                                       2393, 9878, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 25166, 0, 3,
                                                                       22520, 8690, 22736, 2393,
                                                                       2438, 10013, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 25436, 0, 3,
                                                                       22736, 8798, 22952, 2438,
                                                                       2483, 10148, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 25706, 0, 3,
                                                                       22952, 8906, 23168, 2483,
                                                                       2528, 10283, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 25976, 0, 3,
                                                                       23168, 9014, 23384, 2528,
                                                                       2573, 10418, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 26246, 0, 3,
                                                                       23600, 9230, 23816, 2663,
                                                                       2708, 10553, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 26516, 0, 3,
                                                                       23816, 9338, 24032, 2708,
                                                                       2753, 10688, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 26786, 0, 3,
                                                                       24032, 9446, 24248, 2753,
                                                                       2798, 10823, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 27056, 0, 3,
                                                                       24248, 9554, 24464, 2798,
                                                                       2843, 10958, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 27326, 0, 3,
                                                                       24464, 9662, 24680, 2843,
                                                                       2888, 11093, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 27596, 0, 3,
                                                                       24896, 9878, 25166, 2978,
                                                                       3033, 11228, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 27926, 0, 3,
                                                                       25166, 10013, 25436, 3033,
                                                                       3088, 11393, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 28256, 0, 3,
                                                                       25436, 10148, 25706, 3088,
                                                                       3143, 11558, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 28586, 0, 3,
                                                                       25706, 10283, 25976, 3143,
                                                                       3198, 11723, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 28916, 0, 3,
                                                                       26246, 10553, 26516, 3308,
                                                                       3363, 11888, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 29246, 0, 3,
                                                                       26516, 10688, 26786, 3363,
                                                                       3418, 12053, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 29576, 0, 3,
                                                                       26786, 10823, 27056, 3418,
                                                                       3473, 12218, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 29906, 0, 3,
                                                                       27056, 10958, 27326, 3473,
                                                                       3528, 12383, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 30236, 0, 3,
                                                                       27596, 11228, 27926, 3638,
                                                                       3704, 12548, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 30632, 0, 3,
                                                                       27926, 11393, 28256, 3704,
                                                                       3770, 12746, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 31028, 0, 3,
                                                                       28256, 11558, 28586, 3770,
                                                                       3836, 12944, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 31424, 0, 3,
                                                                       28916, 11888, 29246, 3968,
                                                                       4034, 13142, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 31820, 0, 3,
                                                                       29246, 12053, 29576, 4034,
                                                                       4100, 13340, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 32216, 0, 3,
                                                                       29576, 12218, 29906, 4100,
                                                                       4166, 13538, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32612, 3, 4298,
                                                                       4301, 13748, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32622, 3, 4301,
                                                                       4304, 13754, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32632, 3, 4304,
                                                                       4307, 13760, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32642, 3, 4307,
                                                                       4310, 13766, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32652, 3, 4310,
                                                                       4313, 13772, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32662, 3, 4313,
                                                                       4316, 13778, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32672, 3, 4316,
                                                                       4319, 13784, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32682, 3, 4319,
                                                                       4322, 13790, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32692, 3, 4322,
                                                                       4325, 13796, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32702, 3, 4325,
                                                                       4328, 13802, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32712, 3, 4328,
                                                                       4331, 13808, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32722, 3, 4337,
                                                                       4340, 13826, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32732, 3, 4340,
                                                                       4343, 13832, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32742, 3, 4343,
                                                                       4346, 13838, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32752, 3, 4346,
                                                                       4349, 13844, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32762, 3, 4349,
                                                                       4352, 13850, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32772, 3, 4352,
                                                                       4355, 13856, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32782, 3, 4355,
                                                                       4358, 13862, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32792, 3, 4358,
                                                                       4361, 13868, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32802, 3, 4361,
                                                                       4364, 13874, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32812, 3, 4364,
                                                                       4367, 13880, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32822, 3, 4367,
                                                                       4370, 13886, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 32832, 0, 3,
                                                                       32612, 13748, 32622, 4376,
                                                                       4385, 13928, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 32862, 0, 3,
                                                                       32622, 13754, 32632, 4385,
                                                                       4394, 13946, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 32892, 0, 3,
                                                                       32632, 13760, 32642, 4394,
                                                                       4403, 13964, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 32922, 0, 3,
                                                                       32642, 13766, 32652, 4403,
                                                                       4412, 13982, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 32952, 0, 3,
                                                                       32652, 13772, 32662, 4412,
                                                                       4421, 14000, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 32982, 0, 3,
                                                                       32662, 13778, 32672, 4421,
                                                                       4430, 14018, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33012, 0, 3,
                                                                       32672, 13784, 32682, 4430,
                                                                       4439, 14036, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33042, 0, 3,
                                                                       32682, 13790, 32692, 4439,
                                                                       4448, 14054, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33072, 0, 3,
                                                                       32692, 13796, 32702, 4448,
                                                                       4457, 14072, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33102, 0, 3,
                                                                       32702, 13802, 32712, 4457,
                                                                       4466, 14090, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33132, 0, 3,
                                                                       32722, 13826, 32732, 4484,
                                                                       4493, 14144, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33162, 0, 3,
                                                                       32732, 13832, 32742, 4493,
                                                                       4502, 14162, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33192, 0, 3,
                                                                       32742, 13838, 32752, 4502,
                                                                       4511, 14180, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33222, 0, 3,
                                                                       32752, 13844, 32762, 4511,
                                                                       4520, 14198, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33252, 0, 3,
                                                                       32762, 13850, 32772, 4520,
                                                                       4529, 14216, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33282, 0, 3,
                                                                       32772, 13856, 32782, 4529,
                                                                       4538, 14234, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33312, 0, 3,
                                                                       32782, 13862, 32792, 4538,
                                                                       4547, 14252, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33342, 0, 3,
                                                                       32792, 13868, 32802, 4547,
                                                                       4556, 14270, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33372, 0, 3,
                                                                       32802, 13874, 32812, 4556,
                                                                       4565, 14288, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 33402, 0, 3,
                                                                       32812, 13880, 32822, 4565,
                                                                       4574, 14306, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 33432, 0, 3,
                                                                       32832, 13928, 32862, 4592,
                                                                       4610, 14396, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 33492, 0, 3,
                                                                       32862, 13946, 32892, 4610,
                                                                       4628, 14432, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 33552, 0, 3,
                                                                       32892, 13964, 32922, 4628,
                                                                       4646, 14468, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 33612, 0, 3,
                                                                       32922, 13982, 32952, 4646,
                                                                       4664, 14504, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 33672, 0, 3,
                                                                       32952, 14000, 32982, 4664,
                                                                       4682, 14540, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 33732, 0, 3,
                                                                       32982, 14018, 33012, 4682,
                                                                       4700, 14576, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 33792, 0, 3,
                                                                       33012, 14036, 33042, 4700,
                                                                       4718, 14612, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 33852, 0, 3,
                                                                       33042, 14054, 33072, 4718,
                                                                       4736, 14648, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 33912, 0, 3,
                                                                       33072, 14072, 33102, 4736,
                                                                       4754, 14684, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 33972, 0, 3,
                                                                       33132, 14144, 33162, 4790,
                                                                       4808, 14792, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34032, 0, 3,
                                                                       33162, 14162, 33192, 4808,
                                                                       4826, 14828, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34092, 0, 3,
                                                                       33192, 14180, 33222, 4826,
                                                                       4844, 14864, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34152, 0, 3,
                                                                       33222, 14198, 33252, 4844,
                                                                       4862, 14900, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34212, 0, 3,
                                                                       33252, 14216, 33282, 4862,
                                                                       4880, 14936, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34272, 0, 3,
                                                                       33282, 14234, 33312, 4880,
                                                                       4898, 14972, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34332, 0, 3,
                                                                       33312, 14252, 33342, 4898,
                                                                       4916, 15008, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34392, 0, 3,
                                                                       33342, 14270, 33372, 4916,
                                                                       4934, 15044, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 34452, 0, 3,
                                                                       33372, 14288, 33402, 4934,
                                                                       4952, 15080, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 34512, 0, 3,
                                                                       33432, 14396, 33492, 4988,
                                                                       5018, 15236, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 34612, 0, 3,
                                                                       33492, 14432, 33552, 5018,
                                                                       5048, 15296, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 34712, 0, 3,
                                                                       33552, 14468, 33612, 5048,
                                                                       5078, 15356, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 34812, 0, 3,
                                                                       33612, 14504, 33672, 5078,
                                                                       5108, 15416, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 34912, 0, 3,
                                                                       33672, 14540, 33732, 5108,
                                                                       5138, 15476, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 35012, 0, 3,
                                                                       33732, 14576, 33792, 5138,
                                                                       5168, 15536, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 35112, 0, 3,
                                                                       33792, 14612, 33852, 5168,
                                                                       5198, 15596, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 35212, 0, 3,
                                                                       33852, 14648, 33912, 5198,
                                                                       5228, 15656, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 35312, 0, 3,
                                                                       33972, 14792, 34032, 5288,
                                                                       5318, 15836, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 35412, 0, 3,
                                                                       34032, 14828, 34092, 5318,
                                                                       5348, 15896, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 35512, 0, 3,
                                                                       34092, 14864, 34152, 5348,
                                                                       5378, 15956, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 35612, 0, 3,
                                                                       34152, 14900, 34212, 5378,
                                                                       5408, 16016, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 35712, 0, 3,
                                                                       34212, 14936, 34272, 5408,
                                                                       5438, 16076, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 35812, 0, 3,
                                                                       34272, 14972, 34332, 5438,
                                                                       5468, 16136, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 35912, 0, 3,
                                                                       34332, 15008, 34392, 5468,
                                                                       5498, 16196, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 36012, 0, 3,
                                                                       34392, 15044, 34452, 5498,
                                                                       5528, 16256, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 36112, 0, 3,
                                                                       34512, 15236, 34612, 5588,
                                                                       5633, 16496, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 36262, 0, 3,
                                                                       34612, 15296, 34712, 5633,
                                                                       5678, 16586, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 36412, 0, 3,
                                                                       34712, 15356, 34812, 5678,
                                                                       5723, 16676, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 36562, 0, 3,
                                                                       34812, 15416, 34912, 5723,
                                                                       5768, 16766, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 36712, 0, 3,
                                                                       34912, 15476, 35012, 5768,
                                                                       5813, 16856, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 36862, 0, 3,
                                                                       35012, 15536, 35112, 5813,
                                                                       5858, 16946, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 37012, 0, 3,
                                                                       35112, 15596, 35212, 5858,
                                                                       5903, 17036, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 37162, 0, 3,
                                                                       35312, 15836, 35412, 5993,
                                                                       6038, 17306, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 37312, 0, 3,
                                                                       35412, 15896, 35512, 6038,
                                                                       6083, 17396, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 37462, 0, 3,
                                                                       35512, 15956, 35612, 6083,
                                                                       6128, 17486, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 37612, 0, 3,
                                                                       35612, 16016, 35712, 6128,
                                                                       6173, 17576, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 37762, 0, 3,
                                                                       35712, 16076, 35812, 6173,
                                                                       6218, 17666, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 37912, 0, 3,
                                                                       35812, 16136, 35912, 6218,
                                                                       6263, 17756, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 38062, 0, 3,
                                                                       35912, 16196, 36012, 6263,
                                                                       6308, 17846, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 38212, 0, 3,
                                                                       36112, 16496, 36262, 6398,
                                                                       6461, 18188, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 38422, 0, 3,
                                                                       36262, 16586, 36412, 6461,
                                                                       6524, 18314, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 38632, 0, 3,
                                                                       36412, 16676, 36562, 6524,
                                                                       6587, 18440, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 38842, 0, 3,
                                                                       36562, 16766, 36712, 6587,
                                                                       6650, 18566, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 39052, 0, 3,
                                                                       36712, 16856, 36862, 6650,
                                                                       6713, 18692, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 39262, 0, 3,
                                                                       36862, 16946, 37012, 6713,
                                                                       6776, 18818, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 39472, 0, 3,
                                                                       37162, 17306, 37312, 6902,
                                                                       6965, 19196, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 39682, 0, 3,
                                                                       37312, 17396, 37462, 6965,
                                                                       7028, 19322, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 39892, 0, 3,
                                                                       37462, 17486, 37612, 7028,
                                                                       7091, 19448, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 40102, 0, 3,
                                                                       37612, 17576, 37762, 7091,
                                                                       7154, 19574, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 40312, 0, 3,
                                                                       37762, 17666, 37912, 7154,
                                                                       7217, 19700, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 40522, 0, 3,
                                                                       37912, 17756, 38062, 7217,
                                                                       7280, 19826, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 40732, 0, 3,
                                                                       38212, 18188, 38422, 7406,
                                                                       7490, 20288, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 41012, 0, 3,
                                                                       38422, 18314, 38632, 7490,
                                                                       7574, 20456, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 41292, 0, 3,
                                                                       38632, 18440, 38842, 7574,
                                                                       7658, 20624, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 41572, 0, 3,
                                                                       38842, 18566, 39052, 7658,
                                                                       7742, 20792, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 41852, 0, 3,
                                                                       39052, 18692, 39262, 7742,
                                                                       7826, 20960, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 42132, 0, 3,
                                                                       39472, 19196, 39682, 7994,
                                                                       8078, 21464, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 42412, 0, 3,
                                                                       39682, 19322, 39892, 8078,
                                                                       8162, 21632, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 42692, 0, 3,
                                                                       39892, 19448, 40102, 8162,
                                                                       8246, 21800, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 42972, 0, 3,
                                                                       40102, 19574, 40312, 8246,
                                                                       8330, 21968, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 43252, 0, 3,
                                                                       40312, 19700, 40522, 8330,
                                                                       8414, 22136, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 43532, 0, 3,
                                                                       40732, 20288, 41012, 8582,
                                                                       8690, 22736, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 43892, 0, 3,
                                                                       41012, 20456, 41292, 8690,
                                                                       8798, 22952, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 44252, 0, 3,
                                                                       41292, 20624, 41572, 8798,
                                                                       8906, 23168, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 44612, 0, 3,
                                                                       41572, 20792, 41852, 8906,
                                                                       9014, 23384, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 44972, 0, 3,
                                                                       42132, 21464, 42412, 9230,
                                                                       9338, 24032, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 45332, 0, 3,
                                                                       42412, 21632, 42692, 9338,
                                                                       9446, 24248, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 45692, 0, 3,
                                                                       42692, 21800, 42972, 9446,
                                                                       9554, 24464, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 46052, 0, 3,
                                                                       42972, 21968, 43252, 9554,
                                                                       9662, 24680, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 46412, 0, 3,
                                                                       43532, 22736, 43892, 9878,
                                                                       10013, 25436, ncols,
                                                                       gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 46862, 0, 3,
                                                                       43892, 22952, 44252,
                                                                       10013, 10148, 25706,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 47312, 0, 3,
                                                                       44252, 23168, 44612,
                                                                       10148, 10283, 25976,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 47762, 0, 3,
                                                                       44972, 24032, 45332,
                                                                       10553, 10688, 26786,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 48212, 0, 3,
                                                                       45332, 24248, 45692,
                                                                       10688, 10823, 27056,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 48662, 0, 3,
                                                                       45692, 24464, 46052,
                                                                       10823, 10958, 27326,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 49112, 0, 3,
                                                                       46412, 25436, 46862,
                                                                       11228, 11393, 28256,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 49662, 0, 3,
                                                                       46862, 25706, 47312,
                                                                       11393, 11558, 28586,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 50212, 0, 3,
                                                                       47762, 26786, 48212,
                                                                       11888, 12053, 29576,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 50762, 0, 3,
                                                                       48212, 27056, 48662,
                                                                       12053, 12218, 29906,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 51312, 0, 3,
                                                                       49112, 28256, 49662,
                                                                       12548, 12746, 31028,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 51972, 0, 3,
                                                                       50212, 29576, 50762,
                                                                       13142, 13340, 32216,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52632, 3, 13736,
                                                                       13742, 32612, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52647, 3, 13742,
                                                                       13748, 32622, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52662, 3, 13748,
                                                                       13754, 32632, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52677, 3, 13754,
                                                                       13760, 32642, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52692, 3, 13760,
                                                                       13766, 32652, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52707, 3, 13766,
                                                                       13772, 32662, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52722, 3, 13772,
                                                                       13778, 32672, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52737, 3, 13778,
                                                                       13784, 32682, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52752, 3, 13784,
                                                                       13790, 32692, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52767, 3, 13790,
                                                                       13796, 32702, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52782, 3, 13796,
                                                                       13802, 32712, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52797, 3, 13814,
                                                                       13820, 32722, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52812, 3, 13820,
                                                                       13826, 32732, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52827, 3, 13826,
                                                                       13832, 32742, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52842, 3, 13832,
                                                                       13838, 32752, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52857, 3, 13838,
                                                                       13844, 32762, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52872, 3, 13844,
                                                                       13850, 32772, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52887, 3, 13850,
                                                                       13856, 32782, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52902, 3, 13856,
                                                                       13862, 32792, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52917, 3, 13862,
                                                                       13868, 32802, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52932, 3, 13868,
                                                                       13874, 32812, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 52947, 3, 13874,
                                                                       13880, 32822, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 52962, 0, 3,
                                                                       52632, 32612, 52647,
                                                                       13892, 13910, 32832,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53007, 0, 3,
                                                                       52647, 32622, 52662,
                                                                       13910, 13928, 32862,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53052, 0, 3,
                                                                       52662, 32632, 52677,
                                                                       13928, 13946, 32892,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53097, 0, 3,
                                                                       52677, 32642, 52692,
                                                                       13946, 13964, 32922,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53142, 0, 3,
                                                                       52692, 32652, 52707,
                                                                       13964, 13982, 32952,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53187, 0, 3,
                                                                       52707, 32662, 52722,
                                                                       13982, 14000, 32982,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53232, 0, 3,
                                                                       52722, 32672, 52737,
                                                                       14000, 14018, 33012,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53277, 0, 3,
                                                                       52737, 32682, 52752,
                                                                       14018, 14036, 33042,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53322, 0, 3,
                                                                       52752, 32692, 52767,
                                                                       14036, 14054, 33072,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53367, 0, 3,
                                                                       52767, 32702, 52782,
                                                                       14054, 14072, 33102,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53412, 0, 3,
                                                                       52797, 32722, 52812,
                                                                       14108, 14126, 33132,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53457, 0, 3,
                                                                       52812, 32732, 52827,
                                                                       14126, 14144, 33162,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53502, 0, 3,
                                                                       52827, 32742, 52842,
                                                                       14144, 14162, 33192,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53547, 0, 3,
                                                                       52842, 32752, 52857,
                                                                       14162, 14180, 33222,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53592, 0, 3,
                                                                       52857, 32762, 52872,
                                                                       14180, 14198, 33252,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53637, 0, 3,
                                                                       52872, 32772, 52887,
                                                                       14198, 14216, 33282,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53682, 0, 3,
                                                                       52887, 32782, 52902,
                                                                       14216, 14234, 33312,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53727, 0, 3,
                                                                       52902, 32792, 52917,
                                                                       14234, 14252, 33342,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53772, 0, 3,
                                                                       52917, 32802, 52932,
                                                                       14252, 14270, 33372,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 53817, 0, 3,
                                                                       52932, 32812, 52947,
                                                                       14270, 14288, 33402,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 53862, 0, 3,
                                                                       52962, 32832, 53007,
                                                                       14324, 14360, 33432,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 53952, 0, 3,
                                                                       53007, 32862, 53052,
                                                                       14360, 14396, 33492,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 54042, 0, 3,
                                                                       53052, 32892, 53097,
                                                                       14396, 14432, 33552,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 54132, 0, 3,
                                                                       53097, 32922, 53142,
                                                                       14432, 14468, 33612,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 54222, 0, 3,
                                                                       53142, 32952, 53187,
                                                                       14468, 14504, 33672,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 54312, 0, 3,
                                                                       53187, 32982, 53232,
                                                                       14504, 14540, 33732,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 54402, 0, 3,
                                                                       53232, 33012, 53277,
                                                                       14540, 14576, 33792,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 54492, 0, 3,
                                                                       53277, 33042, 53322,
                                                                       14576, 14612, 33852,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 54582, 0, 3,
                                                                       53322, 33072, 53367,
                                                                       14612, 14648, 33912,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 54672, 0, 3,
                                                                       53412, 33132, 53457,
                                                                       14720, 14756, 33972,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 54762, 0, 3,
                                                                       53457, 33162, 53502,
                                                                       14756, 14792, 34032,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 54852, 0, 3,
                                                                       53502, 33192, 53547,
                                                                       14792, 14828, 34092,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 54942, 0, 3,
                                                                       53547, 33222, 53592,
                                                                       14828, 14864, 34152,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 55032, 0, 3,
                                                                       53592, 33252, 53637,
                                                                       14864, 14900, 34212,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 55122, 0, 3,
                                                                       53637, 33282, 53682,
                                                                       14900, 14936, 34272,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 55212, 0, 3,
                                                                       53682, 33312, 53727,
                                                                       14936, 14972, 34332,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 55302, 0, 3,
                                                                       53727, 33342, 53772,
                                                                       14972, 15008, 34392,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 55392, 0, 3,
                                                                       53772, 33372, 53817,
                                                                       15008, 15044, 34452,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 55482, 0, 3,
                                                                       53862, 33432, 53952,
                                                                       15116, 15176, 34512,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 55632, 0, 3,
                                                                       53952, 33492, 54042,
                                                                       15176, 15236, 34612,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 55782, 0, 3,
                                                                       54042, 33552, 54132,
                                                                       15236, 15296, 34712,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 55932, 0, 3,
                                                                       54132, 33612, 54222,
                                                                       15296, 15356, 34812,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 56082, 0, 3,
                                                                       54222, 33672, 54312,
                                                                       15356, 15416, 34912,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 56232, 0, 3,
                                                                       54312, 33732, 54402,
                                                                       15416, 15476, 35012,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 56382, 0, 3,
                                                                       54402, 33792, 54492,
                                                                       15476, 15536, 35112,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 56532, 0, 3,
                                                                       54492, 33852, 54582,
                                                                       15536, 15596, 35212,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 56682, 0, 3,
                                                                       54672, 33972, 54762,
                                                                       15716, 15776, 35312,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 56832, 0, 3,
                                                                       54762, 34032, 54852,
                                                                       15776, 15836, 35412,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 56982, 0, 3,
                                                                       54852, 34092, 54942,
                                                                       15836, 15896, 35512,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 57132, 0, 3,
                                                                       54942, 34152, 55032,
                                                                       15896, 15956, 35612,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 57282, 0, 3,
                                                                       55032, 34212, 55122,
                                                                       15956, 16016, 35712,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 57432, 0, 3,
                                                                       55122, 34272, 55212,
                                                                       16016, 16076, 35812,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 57582, 0, 3,
                                                                       55212, 34332, 55302,
                                                                       16076, 16136, 35912,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 57732, 0, 3,
                                                                       55302, 34392, 55392,
                                                                       16136, 16196, 36012,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 57882, 0, 3,
                                                                       55482, 34512, 55632,
                                                                       16316, 16406, 36112,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 58107, 0, 3,
                                                                       55632, 34612, 55782,
                                                                       16406, 16496, 36262,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 58332, 0, 3,
                                                                       55782, 34712, 55932,
                                                                       16496, 16586, 36412,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 58557, 0, 3,
                                                                       55932, 34812, 56082,
                                                                       16586, 16676, 36562,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 58782, 0, 3,
                                                                       56082, 34912, 56232,
                                                                       16676, 16766, 36712,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 59007, 0, 3,
                                                                       56232, 35012, 56382,
                                                                       16766, 16856, 36862,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 59232, 0, 3,
                                                                       56382, 35112, 56532,
                                                                       16856, 16946, 37012,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 59457, 0, 3,
                                                                       56682, 35312, 56832,
                                                                       17126, 17216, 37162,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 59682, 0, 3,
                                                                       56832, 35412, 56982,
                                                                       17216, 17306, 37312,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 59907, 0, 3,
                                                                       56982, 35512, 57132,
                                                                       17306, 17396, 37462,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 60132, 0, 3,
                                                                       57132, 35612, 57282,
                                                                       17396, 17486, 37612,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 60357, 0, 3,
                                                                       57282, 35712, 57432,
                                                                       17486, 17576, 37762,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 60582, 0, 3,
                                                                       57432, 35812, 57582,
                                                                       17576, 17666, 37912,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 60807, 0, 3,
                                                                       57582, 35912, 57732,
                                                                       17666, 17756, 38062,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 61032, 0, 3,
                                                                       57882, 36112, 58107,
                                                                       17936, 18062, 38212,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 61347, 0, 3,
                                                                       58107, 36262, 58332,
                                                                       18062, 18188, 38422,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 61662, 0, 3,
                                                                       58332, 36412, 58557,
                                                                       18188, 18314, 38632,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 61977, 0, 3,
                                                                       58557, 36562, 58782,
                                                                       18314, 18440, 38842,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 62292, 0, 3,
                                                                       58782, 36712, 59007,
                                                                       18440, 18566, 39052,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 62607, 0, 3,
                                                                       59007, 36862, 59232,
                                                                       18566, 18692, 39262,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 62922, 0, 3,
                                                                       59457, 37162, 59682,
                                                                       18944, 19070, 39472,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 63237, 0, 3,
                                                                       59682, 37312, 59907,
                                                                       19070, 19196, 39682,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 63552, 0, 3,
                                                                       59907, 37462, 60132,
                                                                       19196, 19322, 39892,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 63867, 0, 3,
                                                                       60132, 37612, 60357,
                                                                       19322, 19448, 40102,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 64182, 0, 3,
                                                                       60357, 37762, 60582,
                                                                       19448, 19574, 40312,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 64497, 0, 3,
                                                                       60582, 37912, 60807,
                                                                       19574, 19700, 40522,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 64812, 0, 3,
                                                                       61032, 38212, 61347,
                                                                       19952, 20120, 40732,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 65232, 0, 3,
                                                                       61347, 38422, 61662,
                                                                       20120, 20288, 41012,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 65652, 0, 3,
                                                                       61662, 38632, 61977,
                                                                       20288, 20456, 41292,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 66072, 0, 3,
                                                                       61977, 38842, 62292,
                                                                       20456, 20624, 41572,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 66492, 0, 3,
                                                                       62292, 39052, 62607,
                                                                       20624, 20792, 41852,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 66912, 0, 3,
                                                                       62922, 39472, 63237,
                                                                       21128, 21296, 42132,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 67332, 0, 3,
                                                                       63237, 39682, 63552,
                                                                       21296, 21464, 42412,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 67752, 0, 3,
                                                                       63552, 39892, 63867,
                                                                       21464, 21632, 42692,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 68172, 0, 3,
                                                                       63867, 40102, 64182,
                                                                       21632, 21800, 42972,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 68592, 0, 3,
                                                                       64182, 40312, 64497,
                                                                       21800, 21968, 43252,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 69012, 0, 3,
                                                                       64812, 40732, 65232,
                                                                       22304, 22520, 43532,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 69552, 0, 3,
                                                                       65232, 41012, 65652,
                                                                       22520, 22736, 43892,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 70092, 0, 3,
                                                                       65652, 41292, 66072,
                                                                       22736, 22952, 44252,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 70632, 0, 3,
                                                                       66072, 41572, 66492,
                                                                       22952, 23168, 44612,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 71172, 0, 3,
                                                                       66912, 42132, 67332,
                                                                       23600, 23816, 44972,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 71712, 0, 3,
                                                                       67332, 42412, 67752,
                                                                       23816, 24032, 45332,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 72252, 0, 3,
                                                                       67752, 42692, 68172,
                                                                       24032, 24248, 45692,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 72792, 0, 3,
                                                                       68172, 42972, 68592,
                                                                       24248, 24464, 46052,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 73332, 0, 3,
                                                                       69012, 43532, 69552,
                                                                       24896, 25166, 46412,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 74007, 0, 3,
                                                                       69552, 43892, 70092,
                                                                       25166, 25436, 46862,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 74682, 0, 3,
                                                                       70092, 44252, 70632,
                                                                       25436, 25706, 47312,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 75357, 0, 3,
                                                                       71172, 44972, 71712,
                                                                       26246, 26516, 47762,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 76032, 0, 3,
                                                                       71712, 45332, 72252,
                                                                       26516, 26786, 48212,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 76707, 0, 3,
                                                                       72252, 45692, 72792,
                                                                       26786, 27056, 48662,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 77382, 0, 3,
                                                                       73332, 46412, 74007,
                                                                       27596, 27926, 49112,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 78207, 0, 3,
                                                                       74007, 46862, 74682,
                                                                       27926, 28256, 49662,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 79032, 0, 3,
                                                                       75357, 47762, 76032,
                                                                       28916, 29246, 50212,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 79857, 0, 3,
                                                                       76032, 48212, 76707,
                                                                       29246, 29576, 50762,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 80682, 0, 3,
                                                                       77382, 49112, 78207,
                                                                       30236, 30632, 51312,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 81672, 0, 3,
                                                                       79032, 50212, 79857,
                                                                       31424, 31820, 51972,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 82662, 64812, 420, ncols);

                    simdfunc::contract_primitives(buffer, 83334, 66912, 420, ncols);

                    simdfunc::contract_primitives(buffer, 84006, 69012, 540, ncols);

                    simdfunc::contract_primitives(buffer, 84870, 71172, 540, ncols);

                    simdfunc::contract_primitives(buffer, 85734, 73332, 675, ncols);

                    simdfunc::contract_primitives(buffer, 86814, 75357, 675, ncols);

                    simdfunc::contract_primitives(buffer, 87894, 77382, 825, ncols);

                    simdfunc::contract_primitives(buffer, 89214, 79032, 825, ncols);

                    simdfunc::contract_primitives(buffer, 90534, 80682, 990, ncols);

                    simdfunc::contract_primitives(buffer, 92118, 81672, 990, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 83082, 82662, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 83754, 83334, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 84546, 84006, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 85410, 84870, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 86409, 85734, 45, 1, nmax);

        simdtrf::transform_g_inner(buffer, 87489, 86814, 45, 1, nmax);

        simdtrf::transform_g_inner(buffer, 88719, 87894, 55, 1, nmax);

        simdtrf::transform_g_inner(buffer, 90039, 89214, 55, 1, nmax);

        simdtrf::transform_g_inner(buffer, 91524, 90534, 66, 1, nmax);

        simdtrf::transform_g_inner(buffer, 93108, 92118, 66, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 93702, 83082, 84546, 9, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 94458, 83754, 85410, 9, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 95214, 84546, 86409, 9, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 96186, 85410, 87489, 9, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 97158, 86409, 88719, 9, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 98373, 87489, 90039, 9, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 99588, 88719, 91524, 9, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 101073, 90039, 93108, 9,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 102558, 93702, 95214, 9,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 104070, 94458, 96186, 9,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 105582, 95214, 97158, 9,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 107526, 96186, 98373, 9,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 109470, 97158, 99588, 9,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 111900, 98373, 101073, 9,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 114330, 102558, 105582, 9,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 116850, 104070, 107526, 9,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 119370, 105582, 109470, 9,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 122610, 107526, 111900, 9,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 125850, 114330, 119370, 9,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 129630, 116850, 122610, 9,
                                             nmax);

        simdtrf::transform_g_inner(buffer, 133410, 129630, 28, 9, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 133410, 81, nmax);

        simdtrf::transform_g_inner(buffer, 133410, 125850, 28, 9, nmax);

        simdtrf::transform_i_outer(values + 1053 * nvalues + n * npairs, nvalues, buffer, 133410,
                                   81, nmax);
    }

    for (size_t m = 0; m < 2106; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
