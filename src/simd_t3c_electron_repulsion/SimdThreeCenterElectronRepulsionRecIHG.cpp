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


#include "SimdThreeCenterElectronRepulsionRecIHG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
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
#include "SimdTransferIH.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKF.hpp"
#include "SimdTransferKG.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLD.hpp"
#include "SimdTransferLF.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransferMD.hpp"
#include "SimdTransferMP.hpp"
#include "SimdTransferNP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_ihg_three_center_electron_repulsion(double               *values,
                                            const size_t          npairs,
                                            const size_t          natoms,
                                            const CBasisFunction &a_function,
                                            const CBasisFunction &b_function,
                                            const CBasisFunction &c_function,
                                            const CSimdMatrix    &coordinates,
                                            const CSimdMatrix    &c_coordinates,
                                            CSimdMatrix          &buffer,
                                            const double          threshold) -> void
{
    if (npairs > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ihg_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 104762, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1287 * natoms * npairs, 0.0);

        return;
    }

    const auto pi = mathconst::pi_value();

    // NOTE: a row of the values spans every atom pair of every atom on the ket
    // side, so a kernel handed the block of one atom steps by this to reach the
    // next component -- which is what lets it be the kernel a two-center form
    // uses, unchanged.

    const auto nvalues = natoms * npairs;

    for (size_t n = 0; n < natoms; n++)
    {
        simdfunc::prepare_buffer(buffer, 104762, 55790, 6690, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto mu = a_exps[i] * b_exps[j] / p;

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 15,
                                                             ncols, fj, mu, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 23, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 26, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 29, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 32, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 35, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 38, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 41, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 44, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 47, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 50, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 53, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 56, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 59, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 62, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 65, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 68, 0, 3, 7, 8,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 74, 0, 3, 8, 9,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 80, 0, 3, 9, 10,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 86, 0, 3, 10, 11,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 92, 0, 3, 11, 12,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 98, 0, 3, 12, 13,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 104, 0, 3, 13, 14,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 110, 0, 3, 14, 15,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 116, 0, 3, 15, 16,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 122, 0, 3, 16, 17,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 128, 0, 3, 17, 18,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 134, 0, 3, 18, 19,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 140, 0, 3, 19, 20,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 146, 0, 3, 20, 21,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 152, 0, 3, 23, 26,
                                                                       68, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 162, 0, 3, 26, 29,
                                                                       74, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 172, 0, 3, 29, 32,
                                                                       80, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 182, 0, 3, 32, 35,
                                                                       86, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 192, 0, 3, 35, 38,
                                                                       92, 98, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 202, 0, 3, 38, 41,
                                                                       98, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 212, 0, 3, 41, 44,
                                                                       104, 110, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 222, 0, 3, 44, 47,
                                                                       110, 116, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 232, 0, 3, 47, 50,
                                                                       116, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 242, 0, 3, 50, 53,
                                                                       122, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 252, 0, 3, 53, 56,
                                                                       128, 134, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 262, 0, 3, 56, 59,
                                                                       134, 140, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 272, 0, 3, 59, 62,
                                                                       140, 146, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 282, 0, 3, 68, 74,
                                                                       152, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 297, 0, 3, 74, 80,
                                                                       162, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 312, 0, 3, 80, 86,
                                                                       172, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 327, 0, 3, 86, 92,
                                                                       182, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 342, 0, 3, 92, 98,
                                                                       192, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 357, 0, 3, 98,
                                                                       104, 202, 212, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 372, 0, 3, 104,
                                                                       110, 212, 222, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 387, 0, 3, 110,
                                                                       116, 222, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 402, 0, 3, 116,
                                                                       122, 232, 242, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 417, 0, 3, 122,
                                                                       128, 242, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 432, 0, 3, 128,
                                                                       134, 252, 262, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 447, 0, 3, 134,
                                                                       140, 262, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 462, 0, 3, 152,
                                                                       162, 282, 297, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 483, 0, 3, 162,
                                                                       172, 297, 312, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 504, 0, 3, 172,
                                                                       182, 312, 327, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 525, 0, 3, 182,
                                                                       192, 327, 342, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 546, 0, 3, 192,
                                                                       202, 342, 357, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 567, 0, 3, 202,
                                                                       212, 357, 372, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 588, 0, 3, 212,
                                                                       222, 372, 387, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 609, 0, 3, 222,
                                                                       232, 387, 402, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 630, 0, 3, 232,
                                                                       242, 402, 417, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 651, 0, 3, 242,
                                                                       252, 417, 432, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 672, 0, 3, 252,
                                                                       262, 432, 447, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 693, 0, 3, 282,
                                                                       297, 462, 483, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 721, 0, 3, 297,
                                                                       312, 483, 504, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 749, 0, 3, 312,
                                                                       327, 504, 525, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 777, 0, 3, 327,
                                                                       342, 525, 546, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 805, 0, 3, 342,
                                                                       357, 546, 567, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 833, 0, 3, 357,
                                                                       372, 567, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 861, 0, 3, 372,
                                                                       387, 588, 609, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 889, 0, 3, 387,
                                                                       402, 609, 630, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 917, 0, 3, 402,
                                                                       417, 630, 651, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 945, 0, 3, 417,
                                                                       432, 651, 672, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 973, 0, 3, 462,
                                                                       483, 693, 721, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1009, 0, 3, 483,
                                                                       504, 721, 749, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1045, 0, 3, 504,
                                                                       525, 749, 777, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1081, 0, 3, 525,
                                                                       546, 777, 805, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1117, 0, 3, 546,
                                                                       567, 805, 833, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1153, 0, 3, 567,
                                                                       588, 833, 861, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1189, 0, 3, 588,
                                                                       609, 861, 889, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1225, 0, 3, 609,
                                                                       630, 889, 917, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1261, 0, 3, 630,
                                                                       651, 917, 945, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1297, 0, 3, 693,
                                                                       721, 973, 1009, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1342, 0, 3, 721,
                                                                       749, 1009, 1045, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1387, 0, 3, 749,
                                                                       777, 1045, 1081, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1432, 0, 3, 777,
                                                                       805, 1081, 1117, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1477, 0, 3, 805,
                                                                       833, 1117, 1153, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1522, 0, 3, 833,
                                                                       861, 1153, 1189, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1567, 0, 3, 861,
                                                                       889, 1189, 1225, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1612, 0, 3, 889,
                                                                       917, 1225, 1261, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1657, 0, 3, 973,
                                                                       1009, 1297, 1342, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1712, 0, 3, 1009,
                                                                       1045, 1342, 1387, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1767, 0, 3, 1045,
                                                                       1081, 1387, 1432, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1822, 0, 3, 1081,
                                                                       1117, 1432, 1477, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1877, 0, 3, 1117,
                                                                       1153, 1477, 1522, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1932, 0, 3, 1153,
                                                                       1189, 1522, 1567, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1987, 0, 3, 1189,
                                                                       1225, 1567, 1612, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2042, 0, 3, 1297,
                                                                       1342, 1657, 1712, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2108, 0, 3, 1342,
                                                                       1387, 1712, 1767, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2174, 0, 3, 1387,
                                                                       1432, 1767, 1822, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2240, 0, 3, 1432,
                                                                       1477, 1822, 1877, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2306, 0, 3, 1477,
                                                                       1522, 1877, 1932, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2372, 0, 3, 1522,
                                                                       1567, 1932, 1987, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2438, 0, 3, 1657,
                                                                       1712, 2042, 2108, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2516, 0, 3, 1712,
                                                                       1767, 2108, 2174, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2594, 0, 3, 1767,
                                                                       1822, 2174, 2240, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2672, 0, 3, 1822,
                                                                       1877, 2240, 2306, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2750, 0, 3, 1877,
                                                                       1932, 2306, 2372, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2828, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2831, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2834, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2837, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2840, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2843, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2846, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2849, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2852, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2855, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2858, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2861, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2864, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2867, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2870, 3, 9, 29,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2879, 3, 10, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2888, 3, 11, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2897, 3, 12, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2906, 3, 13, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2915, 3, 14, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2924, 3, 15, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2933, 3, 16, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2942, 3, 17, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2951, 3, 18, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2960, 3, 19, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2969, 3, 20, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2978, 3, 21, 65,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2987, 3, 29, 80,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3005, 3, 32, 86,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3023, 3, 35, 92,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3041, 3, 38, 98,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3059, 3, 41, 104,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3077, 3, 44, 110,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3095, 3, 47, 116,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3113, 3, 50, 122,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3131, 3, 53, 128,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3149, 3, 56, 134,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3167, 3, 59, 140,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3185, 3, 62, 146,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3203, 3, 80, 172,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3233, 3, 86, 182,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3263, 3, 92, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3293, 3, 98, 202,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3323, 3, 104, 212,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3353, 3, 110, 222,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3383, 3, 116, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3413, 3, 122, 242,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3443, 3, 128, 252,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3473, 3, 134, 262,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3503, 3, 140, 272,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3533, 3, 172, 312,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3578, 3, 182, 327,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3623, 3, 192, 342,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3668, 3, 202, 357,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3713, 3, 212, 372,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3758, 3, 222, 387,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3803, 3, 232, 402,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3848, 3, 242, 417,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3893, 3, 252, 432,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3938, 3, 262, 447,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3983, 3, 312, 504,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4046, 3, 327, 525,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4109, 3, 342, 546,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4172, 3, 357, 567,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4235, 3, 372, 588,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4298, 3, 387, 609,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4361, 3, 402, 630,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4424, 3, 417, 651,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4487, 3, 432, 672,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4550, 3, 504, 749,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4634, 3, 525, 777,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4718, 3, 546, 805,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4802, 3, 567, 833,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4886, 3, 588, 861,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4970, 3, 609, 889,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5054, 3, 630, 917,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5138, 3, 651, 945,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5222, 3, 749,
                                                                       1045, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5330, 3, 777,
                                                                       1081, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5438, 3, 805,
                                                                       1117, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5546, 3, 833,
                                                                       1153, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5654, 3, 861,
                                                                       1189, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5762, 3, 889,
                                                                       1225, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5870, 3, 917,
                                                                       1261, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5978, 3, 1045,
                                                                       1387, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6113, 3, 1081,
                                                                       1432, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6248, 3, 1117,
                                                                       1477, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6383, 3, 1153,
                                                                       1522, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6518, 3, 1189,
                                                                       1567, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6653, 3, 1225,
                                                                       1612, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 6788, 3, 1387,
                                                                       1767, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 6953, 3, 1432,
                                                                       1822, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7118, 3, 1477,
                                                                       1877, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7283, 3, 1522,
                                                                       1932, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7448, 3, 1567,
                                                                       1987, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 7613, 3, 1767,
                                                                       2174, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 7811, 3, 1822,
                                                                       2240, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 8009, 3, 1877,
                                                                       2306, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 8207, 3, 1932,
                                                                       2372, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 8405, 3, 2174,
                                                                       2594, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 8639, 3, 2240,
                                                                       2672, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 8873, 3, 2306,
                                                                       2750, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9107, 3, 7, 8,
                                                                       2828, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9113, 3, 8, 9,
                                                                       2831, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9119, 3, 9, 10,
                                                                       2834, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9125, 3, 10, 11,
                                                                       2837, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9131, 3, 11, 12,
                                                                       2840, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9137, 3, 12, 13,
                                                                       2843, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9143, 3, 13, 14,
                                                                       2846, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9149, 3, 14, 15,
                                                                       2849, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9155, 3, 15, 16,
                                                                       2852, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9161, 3, 16, 17,
                                                                       2855, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9167, 3, 17, 18,
                                                                       2858, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9173, 3, 18, 19,
                                                                       2861, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9179, 3, 19, 20,
                                                                       2864, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9185, 3, 20, 21,
                                                                       2867, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9191, 0, 3, 9107,
                                                                       2828, 9113, 2870, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9209, 0, 3, 9113,
                                                                       2831, 9119, 2879, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9227, 0, 3, 9119,
                                                                       2834, 9125, 2888, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9245, 0, 3, 9125,
                                                                       2837, 9131, 2897, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9263, 0, 3, 9131,
                                                                       2840, 9137, 2906, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9281, 0, 3, 9137,
                                                                       2843, 9143, 2915, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9299, 0, 3, 9143,
                                                                       2846, 9149, 2924, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9317, 0, 3, 9149,
                                                                       2849, 9155, 2933, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9335, 0, 3, 9155,
                                                                       2852, 9161, 2942, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9353, 0, 3, 9161,
                                                                       2855, 9167, 2951, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9371, 0, 3, 9167,
                                                                       2858, 9173, 2960, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9389, 0, 3, 9173,
                                                                       2861, 9179, 2969, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9407, 0, 3, 9179,
                                                                       2864, 9185, 2978, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9425, 0, 3, 9191,
                                                                       2870, 9209, 68, 74, 2987,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9461, 0, 3, 9209,
                                                                       2879, 9227, 74, 80, 3005,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9497, 0, 3, 9227,
                                                                       2888, 9245, 80, 86, 3023,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9533, 0, 3, 9245,
                                                                       2897, 9263, 86, 92, 3041,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9569, 0, 3, 9263,
                                                                       2906, 9281, 92, 98, 3059,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9605, 0, 3, 9281,
                                                                       2915, 9299, 98, 104, 3077,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9641, 0, 3, 9299,
                                                                       2924, 9317, 104, 110,
                                                                       3095, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9677, 0, 3, 9317,
                                                                       2933, 9335, 110, 116,
                                                                       3113, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9713, 0, 3, 9335,
                                                                       2942, 9353, 116, 122,
                                                                       3131, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9749, 0, 3, 9353,
                                                                       2951, 9371, 122, 128,
                                                                       3149, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9785, 0, 3, 9371,
                                                                       2960, 9389, 128, 134,
                                                                       3167, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9821, 0, 3, 9389,
                                                                       2969, 9407, 134, 140,
                                                                       3185, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9857, 0, 3, 9425,
                                                                       2987, 9461, 152, 162,
                                                                       3203, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9917, 0, 3, 9461,
                                                                       3005, 9497, 162, 172,
                                                                       3233, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9977, 0, 3, 9497,
                                                                       3023, 9533, 172, 182,
                                                                       3263, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10037, 0, 3, 9533,
                                                                       3041, 9569, 182, 192,
                                                                       3293, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10097, 0, 3, 9569,
                                                                       3059, 9605, 192, 202,
                                                                       3323, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10157, 0, 3, 9605,
                                                                       3077, 9641, 202, 212,
                                                                       3353, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10217, 0, 3, 9641,
                                                                       3095, 9677, 212, 222,
                                                                       3383, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10277, 0, 3, 9677,
                                                                       3113, 9713, 222, 232,
                                                                       3413, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10337, 0, 3, 9713,
                                                                       3131, 9749, 232, 242,
                                                                       3443, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10397, 0, 3, 9749,
                                                                       3149, 9785, 242, 252,
                                                                       3473, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10457, 0, 3, 9785,
                                                                       3167, 9821, 252, 262,
                                                                       3503, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10517, 0, 3, 9857,
                                                                       3203, 9917, 282, 297,
                                                                       3533, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10607, 0, 3, 9917,
                                                                       3233, 9977, 297, 312,
                                                                       3578, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10697, 0, 3, 9977,
                                                                       3263, 10037, 312, 327,
                                                                       3623, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10787, 0, 3,
                                                                       10037, 3293, 10097, 327,
                                                                       342, 3668, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10877, 0, 3,
                                                                       10097, 3323, 10157, 342,
                                                                       357, 3713, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10967, 0, 3,
                                                                       10157, 3353, 10217, 357,
                                                                       372, 3758, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11057, 0, 3,
                                                                       10217, 3383, 10277, 372,
                                                                       387, 3803, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11147, 0, 3,
                                                                       10277, 3413, 10337, 387,
                                                                       402, 3848, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11237, 0, 3,
                                                                       10337, 3443, 10397, 402,
                                                                       417, 3893, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11327, 0, 3,
                                                                       10397, 3473, 10457, 417,
                                                                       432, 3938, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11417, 0, 3,
                                                                       10517, 3533, 10607, 462,
                                                                       483, 3983, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11543, 0, 3,
                                                                       10607, 3578, 10697, 483,
                                                                       504, 4046, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11669, 0, 3,
                                                                       10697, 3623, 10787, 504,
                                                                       525, 4109, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11795, 0, 3,
                                                                       10787, 3668, 10877, 525,
                                                                       546, 4172, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11921, 0, 3,
                                                                       10877, 3713, 10967, 546,
                                                                       567, 4235, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12047, 0, 3,
                                                                       10967, 3758, 11057, 567,
                                                                       588, 4298, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12173, 0, 3,
                                                                       11057, 3803, 11147, 588,
                                                                       609, 4361, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12299, 0, 3,
                                                                       11147, 3848, 11237, 609,
                                                                       630, 4424, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12425, 0, 3,
                                                                       11237, 3893, 11327, 630,
                                                                       651, 4487, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12551, 0, 3,
                                                                       11417, 3983, 11543, 693,
                                                                       721, 4550, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12719, 0, 3,
                                                                       11543, 4046, 11669, 721,
                                                                       749, 4634, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12887, 0, 3,
                                                                       11669, 4109, 11795, 749,
                                                                       777, 4718, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13055, 0, 3,
                                                                       11795, 4172, 11921, 777,
                                                                       805, 4802, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13223, 0, 3,
                                                                       11921, 4235, 12047, 805,
                                                                       833, 4886, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13391, 0, 3,
                                                                       12047, 4298, 12173, 833,
                                                                       861, 4970, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13559, 0, 3,
                                                                       12173, 4361, 12299, 861,
                                                                       889, 5054, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13727, 0, 3,
                                                                       12299, 4424, 12425, 889,
                                                                       917, 5138, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 13895, 0, 3,
                                                                       12551, 4550, 12719, 973,
                                                                       1009, 5222, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 14111, 0, 3,
                                                                       12719, 4634, 12887, 1009,
                                                                       1045, 5330, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 14327, 0, 3,
                                                                       12887, 4718, 13055, 1045,
                                                                       1081, 5438, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 14543, 0, 3,
                                                                       13055, 4802, 13223, 1081,
                                                                       1117, 5546, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 14759, 0, 3,
                                                                       13223, 4886, 13391, 1117,
                                                                       1153, 5654, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 14975, 0, 3,
                                                                       13391, 4970, 13559, 1153,
                                                                       1189, 5762, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15191, 0, 3,
                                                                       13559, 5054, 13727, 1189,
                                                                       1225, 5870, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 15407, 0, 3,
                                                                       13895, 5222, 14111, 1297,
                                                                       1342, 5978, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 15677, 0, 3,
                                                                       14111, 5330, 14327, 1342,
                                                                       1387, 6113, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 15947, 0, 3,
                                                                       14327, 5438, 14543, 1387,
                                                                       1432, 6248, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 16217, 0, 3,
                                                                       14543, 5546, 14759, 1432,
                                                                       1477, 6383, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 16487, 0, 3,
                                                                       14759, 5654, 14975, 1477,
                                                                       1522, 6518, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 16757, 0, 3,
                                                                       14975, 5762, 15191, 1522,
                                                                       1567, 6653, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 17027, 0, 3,
                                                                       15407, 5978, 15677, 1657,
                                                                       1712, 6788, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 17357, 0, 3,
                                                                       15677, 6113, 15947, 1712,
                                                                       1767, 6953, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 17687, 0, 3,
                                                                       15947, 6248, 16217, 1767,
                                                                       1822, 7118, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 18017, 0, 3,
                                                                       16217, 6383, 16487, 1822,
                                                                       1877, 7283, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 18347, 0, 3,
                                                                       16487, 6518, 16757, 1877,
                                                                       1932, 7448, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 18677, 0, 3,
                                                                       17027, 6788, 17357, 2042,
                                                                       2108, 7613, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 19073, 0, 3,
                                                                       17357, 6953, 17687, 2108,
                                                                       2174, 7811, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 19469, 0, 3,
                                                                       17687, 7118, 18017, 2174,
                                                                       2240, 8009, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 19865, 0, 3,
                                                                       18017, 7283, 18347, 2240,
                                                                       2306, 8207, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 20261, 0, 3,
                                                                       18677, 7613, 19073, 2438,
                                                                       2516, 8405, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 20729, 0, 3,
                                                                       19073, 7811, 19469, 2516,
                                                                       2594, 8639, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 21197, 0, 3,
                                                                       19469, 8009, 19865, 2594,
                                                                       2672, 8873, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21665, 3, 2828,
                                                                       2831, 9119, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21675, 3, 2831,
                                                                       2834, 9125, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21685, 3, 2834,
                                                                       2837, 9131, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21695, 3, 2837,
                                                                       2840, 9137, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21705, 3, 2840,
                                                                       2843, 9143, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21715, 3, 2843,
                                                                       2846, 9149, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21725, 3, 2846,
                                                                       2849, 9155, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21735, 3, 2849,
                                                                       2852, 9161, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21745, 3, 2852,
                                                                       2855, 9167, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21755, 3, 2855,
                                                                       2858, 9173, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21765, 3, 2858,
                                                                       2861, 9179, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21775, 3, 2861,
                                                                       2864, 9185, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21785, 0, 3,
                                                                       21665, 9119, 21675, 2870,
                                                                       2879, 9227, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21815, 0, 3,
                                                                       21675, 9125, 21685, 2879,
                                                                       2888, 9245, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21845, 0, 3,
                                                                       21685, 9131, 21695, 2888,
                                                                       2897, 9263, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21875, 0, 3,
                                                                       21695, 9137, 21705, 2897,
                                                                       2906, 9281, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21905, 0, 3,
                                                                       21705, 9143, 21715, 2906,
                                                                       2915, 9299, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21935, 0, 3,
                                                                       21715, 9149, 21725, 2915,
                                                                       2924, 9317, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21965, 0, 3,
                                                                       21725, 9155, 21735, 2924,
                                                                       2933, 9335, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21995, 0, 3,
                                                                       21735, 9161, 21745, 2933,
                                                                       2942, 9353, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22025, 0, 3,
                                                                       21745, 9167, 21755, 2942,
                                                                       2951, 9371, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22055, 0, 3,
                                                                       21755, 9173, 21765, 2951,
                                                                       2960, 9389, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22085, 0, 3,
                                                                       21765, 9179, 21775, 2960,
                                                                       2969, 9407, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22115, 0, 3,
                                                                       21785, 9227, 21815, 2987,
                                                                       3005, 9497, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22175, 0, 3,
                                                                       21815, 9245, 21845, 3005,
                                                                       3023, 9533, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22235, 0, 3,
                                                                       21845, 9263, 21875, 3023,
                                                                       3041, 9569, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22295, 0, 3,
                                                                       21875, 9281, 21905, 3041,
                                                                       3059, 9605, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22355, 0, 3,
                                                                       21905, 9299, 21935, 3059,
                                                                       3077, 9641, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22415, 0, 3,
                                                                       21935, 9317, 21965, 3077,
                                                                       3095, 9677, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22475, 0, 3,
                                                                       21965, 9335, 21995, 3095,
                                                                       3113, 9713, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22535, 0, 3,
                                                                       21995, 9353, 22025, 3113,
                                                                       3131, 9749, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22595, 0, 3,
                                                                       22025, 9371, 22055, 3131,
                                                                       3149, 9785, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22655, 0, 3,
                                                                       22055, 9389, 22085, 3149,
                                                                       3167, 9821, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22715, 0, 3,
                                                                       22115, 9497, 22175, 3203,
                                                                       3233, 9977, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22815, 0, 3,
                                                                       22175, 9533, 22235, 3233,
                                                                       3263, 10037, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22915, 0, 3,
                                                                       22235, 9569, 22295, 3263,
                                                                       3293, 10097, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23015, 0, 3,
                                                                       22295, 9605, 22355, 3293,
                                                                       3323, 10157, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23115, 0, 3,
                                                                       22355, 9641, 22415, 3323,
                                                                       3353, 10217, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23215, 0, 3,
                                                                       22415, 9677, 22475, 3353,
                                                                       3383, 10277, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23315, 0, 3,
                                                                       22475, 9713, 22535, 3383,
                                                                       3413, 10337, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23415, 0, 3,
                                                                       22535, 9749, 22595, 3413,
                                                                       3443, 10397, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23515, 0, 3,
                                                                       22595, 9785, 22655, 3443,
                                                                       3473, 10457, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23615, 0, 3,
                                                                       22715, 9977, 22815, 3533,
                                                                       3578, 10697, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23765, 0, 3,
                                                                       22815, 10037, 22915, 3578,
                                                                       3623, 10787, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23915, 0, 3,
                                                                       22915, 10097, 23015, 3623,
                                                                       3668, 10877, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24065, 0, 3,
                                                                       23015, 10157, 23115, 3668,
                                                                       3713, 10967, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24215, 0, 3,
                                                                       23115, 10217, 23215, 3713,
                                                                       3758, 11057, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24365, 0, 3,
                                                                       23215, 10277, 23315, 3758,
                                                                       3803, 11147, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24515, 0, 3,
                                                                       23315, 10337, 23415, 3803,
                                                                       3848, 11237, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24665, 0, 3,
                                                                       23415, 10397, 23515, 3848,
                                                                       3893, 11327, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 24815, 0, 3,
                                                                       23615, 10697, 23765, 3983,
                                                                       4046, 11669, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25025, 0, 3,
                                                                       23765, 10787, 23915, 4046,
                                                                       4109, 11795, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25235, 0, 3,
                                                                       23915, 10877, 24065, 4109,
                                                                       4172, 11921, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25445, 0, 3,
                                                                       24065, 10967, 24215, 4172,
                                                                       4235, 12047, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25655, 0, 3,
                                                                       24215, 11057, 24365, 4235,
                                                                       4298, 12173, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25865, 0, 3,
                                                                       24365, 11147, 24515, 4298,
                                                                       4361, 12299, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26075, 0, 3,
                                                                       24515, 11237, 24665, 4361,
                                                                       4424, 12425, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 26285, 0, 3,
                                                                       24815, 11669, 25025, 4550,
                                                                       4634, 12887, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 26565, 0, 3,
                                                                       25025, 11795, 25235, 4634,
                                                                       4718, 13055, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 26845, 0, 3,
                                                                       25235, 11921, 25445, 4718,
                                                                       4802, 13223, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27125, 0, 3,
                                                                       25445, 12047, 25655, 4802,
                                                                       4886, 13391, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27405, 0, 3,
                                                                       25655, 12173, 25865, 4886,
                                                                       4970, 13559, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27685, 0, 3,
                                                                       25865, 12299, 26075, 4970,
                                                                       5054, 13727, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 27965, 0, 3,
                                                                       26285, 12887, 26565, 5222,
                                                                       5330, 14327, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 28325, 0, 3,
                                                                       26565, 13055, 26845, 5330,
                                                                       5438, 14543, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 28685, 0, 3,
                                                                       26845, 13223, 27125, 5438,
                                                                       5546, 14759, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 29045, 0, 3,
                                                                       27125, 13391, 27405, 5546,
                                                                       5654, 14975, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 29405, 0, 3,
                                                                       27405, 13559, 27685, 5654,
                                                                       5762, 15191, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 29765, 0, 3,
                                                                       27965, 14327, 28325, 5978,
                                                                       6113, 15947, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 30215, 0, 3,
                                                                       28325, 14543, 28685, 6113,
                                                                       6248, 16217, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 30665, 0, 3,
                                                                       28685, 14759, 29045, 6248,
                                                                       6383, 16487, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 31115, 0, 3,
                                                                       29045, 14975, 29405, 6383,
                                                                       6518, 16757, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 31565, 0, 3,
                                                                       29765, 15947, 30215, 6788,
                                                                       6953, 17687, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 32115, 0, 3,
                                                                       30215, 16217, 30665, 6953,
                                                                       7118, 18017, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 32665, 0, 3,
                                                                       30665, 16487, 31115, 7118,
                                                                       7283, 18347, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 33215, 0, 3,
                                                                       31565, 17687, 32115, 7613,
                                                                       7811, 19469, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 33875, 0, 3,
                                                                       32115, 18017, 32665, 7811,
                                                                       8009, 19865, ncols, gamma,
                                                                       p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 34535, 0, 3,
                                                                       33215, 19469, 33875, 8405,
                                                                       8639, 21197, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35315, 3, 9107,
                                                                       9113, 21665, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35330, 3, 9113,
                                                                       9119, 21675, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35345, 3, 9119,
                                                                       9125, 21685, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35360, 3, 9125,
                                                                       9131, 21695, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35375, 3, 9131,
                                                                       9137, 21705, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35390, 3, 9137,
                                                                       9143, 21715, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35405, 3, 9143,
                                                                       9149, 21725, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35420, 3, 9149,
                                                                       9155, 21735, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35435, 3, 9155,
                                                                       9161, 21745, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35450, 3, 9161,
                                                                       9167, 21755, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35465, 3, 9167,
                                                                       9173, 21765, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35480, 3, 9173,
                                                                       9179, 21775, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35495, 0, 3,
                                                                       35315, 21665, 35330, 9191,
                                                                       9209, 21785, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35540, 0, 3,
                                                                       35330, 21675, 35345, 9209,
                                                                       9227, 21815, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35585, 0, 3,
                                                                       35345, 21685, 35360, 9227,
                                                                       9245, 21845, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35630, 0, 3,
                                                                       35360, 21695, 35375, 9245,
                                                                       9263, 21875, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35675, 0, 3,
                                                                       35375, 21705, 35390, 9263,
                                                                       9281, 21905, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35720, 0, 3,
                                                                       35390, 21715, 35405, 9281,
                                                                       9299, 21935, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35765, 0, 3,
                                                                       35405, 21725, 35420, 9299,
                                                                       9317, 21965, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35810, 0, 3,
                                                                       35420, 21735, 35435, 9317,
                                                                       9335, 21995, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35855, 0, 3,
                                                                       35435, 21745, 35450, 9335,
                                                                       9353, 22025, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35900, 0, 3,
                                                                       35450, 21755, 35465, 9353,
                                                                       9371, 22055, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35945, 0, 3,
                                                                       35465, 21765, 35480, 9371,
                                                                       9389, 22085, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 35990, 0, 3,
                                                                       35495, 21785, 35540, 9425,
                                                                       9461, 22115, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36080, 0, 3,
                                                                       35540, 21815, 35585, 9461,
                                                                       9497, 22175, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36170, 0, 3,
                                                                       35585, 21845, 35630, 9497,
                                                                       9533, 22235, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36260, 0, 3,
                                                                       35630, 21875, 35675, 9533,
                                                                       9569, 22295, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36350, 0, 3,
                                                                       35675, 21905, 35720, 9569,
                                                                       9605, 22355, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36440, 0, 3,
                                                                       35720, 21935, 35765, 9605,
                                                                       9641, 22415, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36530, 0, 3,
                                                                       35765, 21965, 35810, 9641,
                                                                       9677, 22475, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36620, 0, 3,
                                                                       35810, 21995, 35855, 9677,
                                                                       9713, 22535, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36710, 0, 3,
                                                                       35855, 22025, 35900, 9713,
                                                                       9749, 22595, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36800, 0, 3,
                                                                       35900, 22055, 35945, 9749,
                                                                       9785, 22655, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 36890, 0, 3,
                                                                       35990, 22115, 36080, 9857,
                                                                       9917, 22715, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37040, 0, 3,
                                                                       36080, 22175, 36170, 9917,
                                                                       9977, 22815, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37190, 0, 3,
                                                                       36170, 22235, 36260, 9977,
                                                                       10037, 22915, ncols,
                                                                       gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37340, 0, 3,
                                                                       36260, 22295, 36350,
                                                                       10037, 10097, 23015,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37490, 0, 3,
                                                                       36350, 22355, 36440,
                                                                       10097, 10157, 23115,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37640, 0, 3,
                                                                       36440, 22415, 36530,
                                                                       10157, 10217, 23215,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37790, 0, 3,
                                                                       36530, 22475, 36620,
                                                                       10217, 10277, 23315,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37940, 0, 3,
                                                                       36620, 22535, 36710,
                                                                       10277, 10337, 23415,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 38090, 0, 3,
                                                                       36710, 22595, 36800,
                                                                       10337, 10397, 23515,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 38240, 0, 3,
                                                                       36890, 22715, 37040,
                                                                       10517, 10607, 23615,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 38465, 0, 3,
                                                                       37040, 22815, 37190,
                                                                       10607, 10697, 23765,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 38690, 0, 3,
                                                                       37190, 22915, 37340,
                                                                       10697, 10787, 23915,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 38915, 0, 3,
                                                                       37340, 23015, 37490,
                                                                       10787, 10877, 24065,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 39140, 0, 3,
                                                                       37490, 23115, 37640,
                                                                       10877, 10967, 24215,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 39365, 0, 3,
                                                                       37640, 23215, 37790,
                                                                       10967, 11057, 24365,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 39590, 0, 3,
                                                                       37790, 23315, 37940,
                                                                       11057, 11147, 24515,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 39815, 0, 3,
                                                                       37940, 23415, 38090,
                                                                       11147, 11237, 24665,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 40040, 0, 3,
                                                                       38240, 23615, 38465,
                                                                       11417, 11543, 24815,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 40355, 0, 3,
                                                                       38465, 23765, 38690,
                                                                       11543, 11669, 25025,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 40670, 0, 3,
                                                                       38690, 23915, 38915,
                                                                       11669, 11795, 25235,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 40985, 0, 3,
                                                                       38915, 24065, 39140,
                                                                       11795, 11921, 25445,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 41300, 0, 3,
                                                                       39140, 24215, 39365,
                                                                       11921, 12047, 25655,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 41615, 0, 3,
                                                                       39365, 24365, 39590,
                                                                       12047, 12173, 25865,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 41930, 0, 3,
                                                                       39590, 24515, 39815,
                                                                       12173, 12299, 26075,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 42245, 0, 3,
                                                                       40040, 24815, 40355,
                                                                       12551, 12719, 26285,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 42665, 0, 3,
                                                                       40355, 25025, 40670,
                                                                       12719, 12887, 26565,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 43085, 0, 3,
                                                                       40670, 25235, 40985,
                                                                       12887, 13055, 26845,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 43505, 0, 3,
                                                                       40985, 25445, 41300,
                                                                       13055, 13223, 27125,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 43925, 0, 3,
                                                                       41300, 25655, 41615,
                                                                       13223, 13391, 27405,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 44345, 0, 3,
                                                                       41615, 25865, 41930,
                                                                       13391, 13559, 27685,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 44765, 0, 3,
                                                                       42245, 26285, 42665,
                                                                       13895, 14111, 27965,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 45305, 0, 3,
                                                                       42665, 26565, 43085,
                                                                       14111, 14327, 28325,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 45845, 0, 3,
                                                                       43085, 26845, 43505,
                                                                       14327, 14543, 28685,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 46385, 0, 3,
                                                                       43505, 27125, 43925,
                                                                       14543, 14759, 29045,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 46925, 0, 3,
                                                                       43925, 27405, 44345,
                                                                       14759, 14975, 29405,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 47465, 0, 3,
                                                                       44765, 27965, 45305,
                                                                       15407, 15677, 29765,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 48140, 0, 3,
                                                                       45305, 28325, 45845,
                                                                       15677, 15947, 30215,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 48815, 0, 3,
                                                                       45845, 28685, 46385,
                                                                       15947, 16217, 30665,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 49490, 0, 3,
                                                                       46385, 29045, 46925,
                                                                       16217, 16487, 31115,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 50165, 0, 3,
                                                                       47465, 29765, 48140,
                                                                       17027, 17357, 31565,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 50990, 0, 3,
                                                                       48140, 30215, 48815,
                                                                       17357, 17687, 32115,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 51815, 0, 3,
                                                                       48815, 30665, 49490,
                                                                       17687, 18017, 32665,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 52640, 0, 3,
                                                                       50165, 31565, 50990,
                                                                       18677, 19073, 33215,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 53630, 0, 3,
                                                                       50990, 32115, 51815,
                                                                       19073, 19469, 33875,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 54620, 0, 3,
                                                                       52640, 33215, 53630,
                                                                       20261, 20729, 34535,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 55790, 42245, 420, ncols);

                    simdfunc::contract_primitives(buffer, 56462, 44765, 540, ncols);

                    simdfunc::contract_primitives(buffer, 57326, 47465, 675, ncols);

                    simdfunc::contract_primitives(buffer, 58406, 50165, 825, ncols);

                    simdfunc::contract_primitives(buffer, 59726, 52640, 990, ncols);

                    simdfunc::contract_primitives(buffer, 61310, 54620, 1170, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 56210, 55790, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 57002, 56462, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 58001, 57326, 45, 1, nmax);

        simdtrf::transform_g_inner(buffer, 59231, 58406, 55, 1, nmax);

        simdtrf::transform_g_inner(buffer, 60716, 59726, 66, 1, nmax);

        simdtrf::transform_g_inner(buffer, 62480, 61310, 78, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 63182, 56210, 57002, 9, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 63938, 57002, 58001, 9, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 64910, 58001, 59231, 9, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 66125, 59231, 60716, 9, nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 67610, 60716, 62480, 9, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 69392, 63182, 63938, 9, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 70904, 63938, 64910, 9, nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 72848, 64910, 66125, 9, nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 75278, 66125, 67610, 9, nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 78248, 69392, 70904, 9, nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 80768, 70904, 72848, 9, nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 84008, 72848, 75278, 9, nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 88058, 78248, 80768, 9, nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 91838, 80768, 84008, 9, nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 96698, 88058, 91838, 9, nmax);

        simdtrf::transform_h_inner(buffer, 101990, 96698, 28, 9, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 101990, 99, nmax);
    }

    for (size_t m = 0; m < 1287; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
