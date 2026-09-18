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


#include "SimdThreeCenterElectronRepulsionRsRecIID.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
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
#include "SimdTransformD.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_iid_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_iid_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 125896, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1690 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 125896, 38228, 8323, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5468, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5471, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5474, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5477, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5480, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5483, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5486, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5489, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5492, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5495, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5498, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5501, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5504, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5507, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5510, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5513, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5516, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5519, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5522, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5525, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5528, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5531, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5534, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5537, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5540, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5543, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5546, 3, 9, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5555, 3, 10, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5564, 3, 11, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5573, 3, 12, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5582, 3, 13, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5591, 3, 14, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5600, 3, 15, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5609, 3, 16, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5618, 3, 17, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5627, 3, 18, 71,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5636, 3, 19, 74,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5645, 3, 20, 77,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5654, 3, 25, 86,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5663, 3, 26, 89,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5672, 3, 27, 92,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5681, 3, 28, 95,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5690, 3, 29, 98,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5699, 3, 30, 101,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5708, 3, 31, 104,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5717, 3, 32, 107,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5726, 3, 33, 110,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5735, 3, 34, 113,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5744, 3, 35, 116,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5753, 3, 36, 119,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5762, 3, 44, 134,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5780, 3, 47, 140,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5798, 3, 50, 146,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5816, 3, 53, 152,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5834, 3, 56, 158,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5852, 3, 59, 164,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5870, 3, 62, 170,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5888, 3, 65, 176,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5906, 3, 68, 182,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5924, 3, 71, 188,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5942, 3, 74, 194,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5960, 3, 86, 212,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5978, 3, 89, 218,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5996, 3, 92, 224,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6014, 3, 95, 230,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6032, 3, 98, 236,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6050, 3, 101, 242,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6068, 3, 104, 248,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6086, 3, 107, 254,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6104, 3, 110, 260,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6122, 3, 113, 266,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 6140, 3, 116, 272,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6158, 3, 134, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6188, 3, 140, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6218, 3, 146, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6248, 3, 152, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6278, 3, 158, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6308, 3, 164, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6338, 3, 170, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6368, 3, 176, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6398, 3, 182, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6428, 3, 188, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6458, 3, 212, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6488, 3, 218, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6518, 3, 224, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6548, 3, 230, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6578, 3, 236, 458,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6608, 3, 242, 468,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6638, 3, 248, 478,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6668, 3, 254, 488,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6698, 3, 260, 498,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6728, 3, 266, 508,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6758, 3, 298, 548,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6803, 3, 308, 563,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6848, 3, 318, 578,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6893, 3, 328, 593,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6938, 3, 338, 608,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6983, 3, 348, 623,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7028, 3, 358, 638,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7073, 3, 368, 653,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7118, 3, 378, 668,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7163, 3, 418, 713,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7208, 3, 428, 728,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7253, 3, 438, 743,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7298, 3, 448, 758,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7343, 3, 458, 773,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7388, 3, 468, 788,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7433, 3, 478, 803,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7478, 3, 488, 818,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7523, 3, 498, 833,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7568, 3, 548, 890,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7631, 3, 563, 911,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7694, 3, 578, 932,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7757, 3, 593, 953,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7820, 3, 608, 974,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7883, 3, 623, 995,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7946, 3, 638,
                                                                       1016, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8009, 3, 653,
                                                                       1037, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8072, 3, 713,
                                                                       1100, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8135, 3, 728,
                                                                       1121, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8198, 3, 743,
                                                                       1142, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8261, 3, 758,
                                                                       1163, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8324, 3, 773,
                                                                       1184, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8387, 3, 788,
                                                                       1205, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8450, 3, 803,
                                                                       1226, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 8513, 3, 818,
                                                                       1247, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8576, 3, 890,
                                                                       1324, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8660, 3, 911,
                                                                       1352, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8744, 3, 932,
                                                                       1380, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8828, 3, 953,
                                                                       1408, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8912, 3, 974,
                                                                       1436, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8996, 3, 995,
                                                                       1464, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9080, 3, 1016,
                                                                       1492, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9164, 3, 1100,
                                                                       1576, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9248, 3, 1121,
                                                                       1604, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9332, 3, 1142,
                                                                       1632, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9416, 3, 1163,
                                                                       1660, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9500, 3, 1184,
                                                                       1688, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9584, 3, 1205,
                                                                       1716, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9668, 3, 1226,
                                                                       1744, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9752, 3, 1324,
                                                                       1844, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9860, 3, 1352,
                                                                       1880, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9968, 3, 1380,
                                                                       1916, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10076, 3, 1408,
                                                                       1952, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10184, 3, 1436,
                                                                       1988, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10292, 3, 1464,
                                                                       2024, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10400, 3, 1576,
                                                                       2132, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10508, 3, 1604,
                                                                       2168, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10616, 3, 1632,
                                                                       2204, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10724, 3, 1660,
                                                                       2240, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10832, 3, 1688,
                                                                       2276, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10940, 3, 1716,
                                                                       2312, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11048, 3, 1844,
                                                                       2438, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11183, 3, 1880,
                                                                       2483, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11318, 3, 1916,
                                                                       2528, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11453, 3, 1952,
                                                                       2573, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11588, 3, 1988,
                                                                       2618, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11723, 3, 2132,
                                                                       2753, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11858, 3, 2168,
                                                                       2798, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11993, 3, 2204,
                                                                       2843, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12128, 3, 2240,
                                                                       2888, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 12263, 3, 2276,
                                                                       2933, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12398, 3, 2438,
                                                                       3088, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12563, 3, 2483,
                                                                       3143, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12728, 3, 2528,
                                                                       3198, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12893, 3, 2573,
                                                                       3253, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13058, 3, 2753,
                                                                       3418, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13223, 3, 2798,
                                                                       3473, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13388, 3, 2843,
                                                                       3528, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13553, 3, 2888,
                                                                       3583, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13718, 3, 3088,
                                                                       3770, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13916, 3, 3143,
                                                                       3836, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14114, 3, 3198,
                                                                       3902, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14312, 3, 3418,
                                                                       4100, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14510, 3, 3473,
                                                                       4166, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14708, 3, 3528,
                                                                       4232, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 14906, 3, 3770,
                                                                       4454, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 15140, 3, 3836,
                                                                       4532, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 15374, 3, 4100,
                                                                       4766, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 15608, 3, 4166,
                                                                       4844, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 15842, 3, 4454,
                                                                       5104, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 16115, 3, 4766,
                                                                       5377, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16388, 3, 7, 8,
                                                                       5468, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16394, 3, 8, 9,
                                                                       5471, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16400, 3, 9, 10,
                                                                       5474, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16406, 3, 10, 11,
                                                                       5477, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16412, 3, 11, 12,
                                                                       5480, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16418, 3, 12, 13,
                                                                       5483, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16424, 3, 13, 14,
                                                                       5486, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16430, 3, 14, 15,
                                                                       5489, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16436, 3, 15, 16,
                                                                       5492, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16442, 3, 16, 17,
                                                                       5495, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16448, 3, 17, 18,
                                                                       5498, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16454, 3, 18, 19,
                                                                       5501, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16460, 3, 19, 20,
                                                                       5504, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16466, 3, 23, 24,
                                                                       5507, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16472, 3, 24, 25,
                                                                       5510, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16478, 3, 25, 26,
                                                                       5513, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16484, 3, 26, 27,
                                                                       5516, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16490, 3, 27, 28,
                                                                       5519, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16496, 3, 28, 29,
                                                                       5522, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16502, 3, 29, 30,
                                                                       5525, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16508, 3, 30, 31,
                                                                       5528, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16514, 3, 31, 32,
                                                                       5531, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16520, 3, 32, 33,
                                                                       5534, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16526, 3, 33, 34,
                                                                       5537, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16532, 3, 34, 35,
                                                                       5540, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16538, 3, 35, 36,
                                                                       5543, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16544, 0, 3,
                                                                       16388, 5468, 16394, 5546,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16562, 0, 3,
                                                                       16394, 5471, 16400, 5555,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16580, 0, 3,
                                                                       16400, 5474, 16406, 5564,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16598, 0, 3,
                                                                       16406, 5477, 16412, 5573,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16616, 0, 3,
                                                                       16412, 5480, 16418, 5582,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16634, 0, 3,
                                                                       16418, 5483, 16424, 5591,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16652, 0, 3,
                                                                       16424, 5486, 16430, 5600,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16670, 0, 3,
                                                                       16430, 5489, 16436, 5609,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16688, 0, 3,
                                                                       16436, 5492, 16442, 5618,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16706, 0, 3,
                                                                       16442, 5495, 16448, 5627,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16724, 0, 3,
                                                                       16448, 5498, 16454, 5636,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16742, 0, 3,
                                                                       16454, 5501, 16460, 5645,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16760, 0, 3,
                                                                       16466, 5507, 16472, 5654,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16778, 0, 3,
                                                                       16472, 5510, 16478, 5663,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16796, 0, 3,
                                                                       16478, 5513, 16484, 5672,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16814, 0, 3,
                                                                       16484, 5516, 16490, 5681,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16832, 0, 3,
                                                                       16490, 5519, 16496, 5690,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16850, 0, 3,
                                                                       16496, 5522, 16502, 5699,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16868, 0, 3,
                                                                       16502, 5525, 16508, 5708,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16886, 0, 3,
                                                                       16508, 5528, 16514, 5717,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16904, 0, 3,
                                                                       16514, 5531, 16520, 5726,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16922, 0, 3,
                                                                       16520, 5534, 16526, 5735,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16940, 0, 3,
                                                                       16526, 5537, 16532, 5744,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 16958, 0, 3,
                                                                       16532, 5540, 16538, 5753,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16976, 0, 3,
                                                                       16544, 5546, 16562, 122,
                                                                       128, 5762, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17012, 0, 3,
                                                                       16562, 5555, 16580, 128,
                                                                       134, 5780, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17048, 0, 3,
                                                                       16580, 5564, 16598, 134,
                                                                       140, 5798, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17084, 0, 3,
                                                                       16598, 5573, 16616, 140,
                                                                       146, 5816, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17120, 0, 3,
                                                                       16616, 5582, 16634, 146,
                                                                       152, 5834, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17156, 0, 3,
                                                                       16634, 5591, 16652, 152,
                                                                       158, 5852, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17192, 0, 3,
                                                                       16652, 5600, 16670, 158,
                                                                       164, 5870, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17228, 0, 3,
                                                                       16670, 5609, 16688, 164,
                                                                       170, 5888, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17264, 0, 3,
                                                                       16688, 5618, 16706, 170,
                                                                       176, 5906, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17300, 0, 3,
                                                                       16706, 5627, 16724, 176,
                                                                       182, 5924, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17336, 0, 3,
                                                                       16724, 5636, 16742, 182,
                                                                       188, 5942, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17372, 0, 3,
                                                                       16760, 5654, 16778, 200,
                                                                       206, 5960, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17408, 0, 3,
                                                                       16778, 5663, 16796, 206,
                                                                       212, 5978, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17444, 0, 3,
                                                                       16796, 5672, 16814, 212,
                                                                       218, 5996, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17480, 0, 3,
                                                                       16814, 5681, 16832, 218,
                                                                       224, 6014, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17516, 0, 3,
                                                                       16832, 5690, 16850, 224,
                                                                       230, 6032, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17552, 0, 3,
                                                                       16850, 5699, 16868, 230,
                                                                       236, 6050, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17588, 0, 3,
                                                                       16868, 5708, 16886, 236,
                                                                       242, 6068, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17624, 0, 3,
                                                                       16886, 5717, 16904, 242,
                                                                       248, 6086, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17660, 0, 3,
                                                                       16904, 5726, 16922, 248,
                                                                       254, 6104, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17696, 0, 3,
                                                                       16922, 5735, 16940, 254,
                                                                       260, 6122, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 17732, 0, 3,
                                                                       16940, 5744, 16958, 260,
                                                                       266, 6140, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 17768, 0, 3,
                                                                       16976, 5762, 17012, 278,
                                                                       288, 6158, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 17828, 0, 3,
                                                                       17012, 5780, 17048, 288,
                                                                       298, 6188, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 17888, 0, 3,
                                                                       17048, 5798, 17084, 298,
                                                                       308, 6218, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 17948, 0, 3,
                                                                       17084, 5816, 17120, 308,
                                                                       318, 6248, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18008, 0, 3,
                                                                       17120, 5834, 17156, 318,
                                                                       328, 6278, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18068, 0, 3,
                                                                       17156, 5852, 17192, 328,
                                                                       338, 6308, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18128, 0, 3,
                                                                       17192, 5870, 17228, 338,
                                                                       348, 6338, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18188, 0, 3,
                                                                       17228, 5888, 17264, 348,
                                                                       358, 6368, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18248, 0, 3,
                                                                       17264, 5906, 17300, 358,
                                                                       368, 6398, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18308, 0, 3,
                                                                       17300, 5924, 17336, 368,
                                                                       378, 6428, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18368, 0, 3,
                                                                       17372, 5960, 17408, 398,
                                                                       408, 6458, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18428, 0, 3,
                                                                       17408, 5978, 17444, 408,
                                                                       418, 6488, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18488, 0, 3,
                                                                       17444, 5996, 17480, 418,
                                                                       428, 6518, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18548, 0, 3,
                                                                       17480, 6014, 17516, 428,
                                                                       438, 6548, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18608, 0, 3,
                                                                       17516, 6032, 17552, 438,
                                                                       448, 6578, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18668, 0, 3,
                                                                       17552, 6050, 17588, 448,
                                                                       458, 6608, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18728, 0, 3,
                                                                       17588, 6068, 17624, 458,
                                                                       468, 6638, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18788, 0, 3,
                                                                       17624, 6086, 17660, 468,
                                                                       478, 6668, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18848, 0, 3,
                                                                       17660, 6104, 17696, 478,
                                                                       488, 6698, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 18908, 0, 3,
                                                                       17696, 6122, 17732, 488,
                                                                       498, 6728, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 18968, 0, 3,
                                                                       17768, 6158, 17828, 518,
                                                                       533, 6758, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19058, 0, 3,
                                                                       17828, 6188, 17888, 533,
                                                                       548, 6803, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19148, 0, 3,
                                                                       17888, 6218, 17948, 548,
                                                                       563, 6848, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19238, 0, 3,
                                                                       17948, 6248, 18008, 563,
                                                                       578, 6893, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19328, 0, 3,
                                                                       18008, 6278, 18068, 578,
                                                                       593, 6938, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19418, 0, 3,
                                                                       18068, 6308, 18128, 593,
                                                                       608, 6983, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19508, 0, 3,
                                                                       18128, 6338, 18188, 608,
                                                                       623, 7028, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19598, 0, 3,
                                                                       18188, 6368, 18248, 623,
                                                                       638, 7073, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19688, 0, 3,
                                                                       18248, 6398, 18308, 638,
                                                                       653, 7118, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19778, 0, 3,
                                                                       18368, 6458, 18428, 683,
                                                                       698, 7163, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19868, 0, 3,
                                                                       18428, 6488, 18488, 698,
                                                                       713, 7208, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19958, 0, 3,
                                                                       18488, 6518, 18548, 713,
                                                                       728, 7253, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20048, 0, 3,
                                                                       18548, 6548, 18608, 728,
                                                                       743, 7298, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20138, 0, 3,
                                                                       18608, 6578, 18668, 743,
                                                                       758, 7343, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20228, 0, 3,
                                                                       18668, 6608, 18728, 758,
                                                                       773, 7388, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20318, 0, 3,
                                                                       18728, 6638, 18788, 773,
                                                                       788, 7433, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20408, 0, 3,
                                                                       18788, 6668, 18848, 788,
                                                                       803, 7478, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20498, 0, 3,
                                                                       18848, 6698, 18908, 803,
                                                                       818, 7523, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 20588, 0, 3,
                                                                       18968, 6758, 19058, 848,
                                                                       869, 7568, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 20714, 0, 3,
                                                                       19058, 6803, 19148, 869,
                                                                       890, 7631, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 20840, 0, 3,
                                                                       19148, 6848, 19238, 890,
                                                                       911, 7694, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 20966, 0, 3,
                                                                       19238, 6893, 19328, 911,
                                                                       932, 7757, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21092, 0, 3,
                                                                       19328, 6938, 19418, 932,
                                                                       953, 7820, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21218, 0, 3,
                                                                       19418, 6983, 19508, 953,
                                                                       974, 7883, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21344, 0, 3,
                                                                       19508, 7028, 19598, 974,
                                                                       995, 7946, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21470, 0, 3,
                                                                       19598, 7073, 19688, 995,
                                                                       1016, 8009, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21596, 0, 3,
                                                                       19778, 7163, 19868, 1058,
                                                                       1079, 8072, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21722, 0, 3,
                                                                       19868, 7208, 19958, 1079,
                                                                       1100, 8135, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21848, 0, 3,
                                                                       19958, 7253, 20048, 1100,
                                                                       1121, 8198, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21974, 0, 3,
                                                                       20048, 7298, 20138, 1121,
                                                                       1142, 8261, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22100, 0, 3,
                                                                       20138, 7343, 20228, 1142,
                                                                       1163, 8324, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22226, 0, 3,
                                                                       20228, 7388, 20318, 1163,
                                                                       1184, 8387, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22352, 0, 3,
                                                                       20318, 7433, 20408, 1184,
                                                                       1205, 8450, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22478, 0, 3,
                                                                       20408, 7478, 20498, 1205,
                                                                       1226, 8513, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22604, 0, 3,
                                                                       20588, 7568, 20714, 1268,
                                                                       1296, 8576, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22772, 0, 3,
                                                                       20714, 7631, 20840, 1296,
                                                                       1324, 8660, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22940, 0, 3,
                                                                       20840, 7694, 20966, 1324,
                                                                       1352, 8744, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23108, 0, 3,
                                                                       20966, 7757, 21092, 1352,
                                                                       1380, 8828, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23276, 0, 3,
                                                                       21092, 7820, 21218, 1380,
                                                                       1408, 8912, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23444, 0, 3,
                                                                       21218, 7883, 21344, 1408,
                                                                       1436, 8996, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23612, 0, 3,
                                                                       21344, 7946, 21470, 1436,
                                                                       1464, 9080, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23780, 0, 3,
                                                                       21596, 8072, 21722, 1520,
                                                                       1548, 9164, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23948, 0, 3,
                                                                       21722, 8135, 21848, 1548,
                                                                       1576, 9248, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24116, 0, 3,
                                                                       21848, 8198, 21974, 1576,
                                                                       1604, 9332, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24284, 0, 3,
                                                                       21974, 8261, 22100, 1604,
                                                                       1632, 9416, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24452, 0, 3,
                                                                       22100, 8324, 22226, 1632,
                                                                       1660, 9500, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24620, 0, 3,
                                                                       22226, 8387, 22352, 1660,
                                                                       1688, 9584, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24788, 0, 3,
                                                                       22352, 8450, 22478, 1688,
                                                                       1716, 9668, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 24956, 0, 3,
                                                                       22604, 8576, 22772, 1772,
                                                                       1808, 9752, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25172, 0, 3,
                                                                       22772, 8660, 22940, 1808,
                                                                       1844, 9860, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25388, 0, 3,
                                                                       22940, 8744, 23108, 1844,
                                                                       1880, 9968, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25604, 0, 3,
                                                                       23108, 8828, 23276, 1880,
                                                                       1916, 10076, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25820, 0, 3,
                                                                       23276, 8912, 23444, 1916,
                                                                       1952, 10184, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26036, 0, 3,
                                                                       23444, 8996, 23612, 1952,
                                                                       1988, 10292, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26252, 0, 3,
                                                                       23780, 9164, 23948, 2060,
                                                                       2096, 10400, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26468, 0, 3,
                                                                       23948, 9248, 24116, 2096,
                                                                       2132, 10508, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26684, 0, 3,
                                                                       24116, 9332, 24284, 2132,
                                                                       2168, 10616, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26900, 0, 3,
                                                                       24284, 9416, 24452, 2168,
                                                                       2204, 10724, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27116, 0, 3,
                                                                       24452, 9500, 24620, 2204,
                                                                       2240, 10832, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27332, 0, 3,
                                                                       24620, 9584, 24788, 2240,
                                                                       2276, 10940, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 27548, 0, 3,
                                                                       24956, 9752, 25172, 2348,
                                                                       2393, 11048, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 27818, 0, 3,
                                                                       25172, 9860, 25388, 2393,
                                                                       2438, 11183, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28088, 0, 3,
                                                                       25388, 9968, 25604, 2438,
                                                                       2483, 11318, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28358, 0, 3,
                                                                       25604, 10076, 25820, 2483,
                                                                       2528, 11453, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28628, 0, 3,
                                                                       25820, 10184, 26036, 2528,
                                                                       2573, 11588, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28898, 0, 3,
                                                                       26252, 10400, 26468, 2663,
                                                                       2708, 11723, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29168, 0, 3,
                                                                       26468, 10508, 26684, 2708,
                                                                       2753, 11858, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29438, 0, 3,
                                                                       26684, 10616, 26900, 2753,
                                                                       2798, 11993, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29708, 0, 3,
                                                                       26900, 10724, 27116, 2798,
                                                                       2843, 12128, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29978, 0, 3,
                                                                       27116, 10832, 27332, 2843,
                                                                       2888, 12263, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 30248, 0, 3,
                                                                       27548, 11048, 27818, 2978,
                                                                       3033, 12398, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 30578, 0, 3,
                                                                       27818, 11183, 28088, 3033,
                                                                       3088, 12563, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 30908, 0, 3,
                                                                       28088, 11318, 28358, 3088,
                                                                       3143, 12728, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 31238, 0, 3,
                                                                       28358, 11453, 28628, 3143,
                                                                       3198, 12893, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 31568, 0, 3,
                                                                       28898, 11723, 29168, 3308,
                                                                       3363, 13058, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 31898, 0, 3,
                                                                       29168, 11858, 29438, 3363,
                                                                       3418, 13223, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 32228, 0, 3,
                                                                       29438, 11993, 29708, 3418,
                                                                       3473, 13388, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 32558, 0, 3,
                                                                       29708, 12128, 29978, 3473,
                                                                       3528, 13553, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 32888, 0, 3,
                                                                       30248, 12398, 30578, 3638,
                                                                       3704, 13718, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 33284, 0, 3,
                                                                       30578, 12563, 30908, 3704,
                                                                       3770, 13916, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 33680, 0, 3,
                                                                       30908, 12728, 31238, 3770,
                                                                       3836, 14114, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 34076, 0, 3,
                                                                       31568, 13058, 31898, 3968,
                                                                       4034, 14312, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 34472, 0, 3,
                                                                       31898, 13223, 32228, 4034,
                                                                       4100, 14510, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 34868, 0, 3,
                                                                       32228, 13388, 32558, 4100,
                                                                       4166, 14708, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 35264, 0, 3,
                                                                       32888, 13718, 33284, 4298,
                                                                       4376, 14906, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 35732, 0, 3,
                                                                       33284, 13916, 33680, 4376,
                                                                       4454, 15140, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 36200, 0, 3,
                                                                       34076, 14312, 34472, 4610,
                                                                       4688, 15374, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 36668, 0, 3,
                                                                       34472, 14510, 34868, 4688,
                                                                       4766, 15608, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 37136, 0, 3,
                                                                       35264, 14906, 35732, 4922,
                                                                       5013, 15842, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 37682, 0, 3,
                                                                       36200, 15374, 36668, 5195,
                                                                       5286, 16115, ncols, gamma,
                                                                       p, q);

                    simdfunc::contract_primitives(buffer, 38228, 22604, 168, ncols);

                    simdfunc::contract_primitives(buffer, 38536, 23780, 168, ncols);

                    simdfunc::contract_primitives(buffer, 38844, 24956, 216, ncols);

                    simdfunc::contract_primitives(buffer, 39240, 26252, 216, ncols);

                    simdfunc::contract_primitives(buffer, 39636, 27548, 270, ncols);

                    simdfunc::contract_primitives(buffer, 40131, 28898, 270, ncols);

                    simdfunc::contract_primitives(buffer, 40626, 30248, 330, ncols);

                    simdfunc::contract_primitives(buffer, 41231, 31568, 330, ncols);

                    simdfunc::contract_primitives(buffer, 41836, 32888, 396, ncols);

                    simdfunc::contract_primitives(buffer, 42562, 34076, 396, ncols);

                    simdfunc::contract_primitives(buffer, 43288, 35264, 468, ncols);

                    simdfunc::contract_primitives(buffer, 44146, 36200, 468, ncols);

                    simdfunc::contract_primitives(buffer, 45004, 37136, 546, ncols);

                    simdfunc::contract_primitives(buffer, 46005, 37682, 546, ncols);
                }
            }
        }

        simdtrf::transform_d_inner(buffer, 38396, 38228, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 38704, 38536, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 39060, 38844, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 39456, 39240, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 39906, 39636, 45, 1, nmax);

        simdtrf::transform_d_inner(buffer, 40401, 40131, 45, 1, nmax);

        simdtrf::transform_d_inner(buffer, 40956, 40626, 55, 1, nmax);

        simdtrf::transform_d_inner(buffer, 41561, 41231, 55, 1, nmax);

        simdtrf::transform_d_inner(buffer, 42232, 41836, 66, 1, nmax);

        simdtrf::transform_d_inner(buffer, 42958, 42562, 66, 1, nmax);

        simdtrf::transform_d_inner(buffer, 43756, 43288, 78, 1, nmax);

        simdtrf::transform_d_inner(buffer, 44614, 44146, 78, 1, nmax);

        simdtrf::transform_d_inner(buffer, 45550, 45004, 91, 1, nmax);

        simdtrf::transform_d_inner(buffer, 46551, 46005, 91, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 47006, 38396, 39060, 5, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 47426, 38704, 39456, 5, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 47846, 39060, 39906, 5, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 48386, 39456, 40401, 5, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 48926, 39906, 40956, 5, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 49601, 40401, 41561, 5, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 50276, 40956, 42232, 5, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 51101, 41561, 42958, 5, nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 51926, 42232, 43756, 5, nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 52916, 42958, 44614, 5, nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 53906, 43756, 45550, 5, nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 55076, 44614, 46551, 5, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 56246, 47006, 47846, 5, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 57086, 47426, 48386, 5, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 57926, 47846, 48926, 5, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 59006, 48386, 49601, 5, nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 60086, 48926, 50276, 5, nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 61436, 49601, 51101, 5, nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 62786, 50276, 51926, 5, nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 64436, 51101, 52916, 5, nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 66086, 51926, 53906, 5, nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 68066, 52916, 55076, 5, nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 70046, 56246, 57926, 5, nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 71446, 57086, 59006, 5, nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 72846, 57926, 60086, 5, nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 74646, 59006, 61436, 5, nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 76446, 60086, 62786, 5, nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 78696, 61436, 64436, 5, nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 80946, 62786, 66086, 5, nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 83696, 64436, 68066, 5, nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 86446, 70046, 72846, 5, nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 88546, 71446, 74646, 5, nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 90646, 72846, 76446, 5, nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 93346, 74646, 78696, 5, nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 96046, 76446, 80946, 5, nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 99421, 78696, 83696, 5, nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 102796, 86446, 90646, 5,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 105736, 88546, 93346, 5,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 108676, 90646, 96046, 5,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 112456, 93346, 99421, 5,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 116236, 102796, 108676, 5,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 120156, 105736, 112456, 5,
                                             nmax);

        simdtrf::transform_i_inner(buffer, 124076, 120156, 28, 5, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 124076, 65, nmax);

        simdtrf::transform_i_inner(buffer, 124076, 116236, 28, 5, nmax);

        simdtrf::transform_i_outer(values + 845 * nvalues + n * npairs, nvalues, buffer, 124076,
                                   65, nmax);
    }

    for (size_t m = 0; m < 1690; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
