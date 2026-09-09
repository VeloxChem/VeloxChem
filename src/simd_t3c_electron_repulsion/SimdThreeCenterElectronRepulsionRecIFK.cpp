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


#include "SimdThreeCenterElectronRepulsionRecIFK.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
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
#include "SimdTransferID.hpp"
#include "SimdTransferIF.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_ifk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_ifk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 147176, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1365 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 147176, 121007, 7539, dimensions);

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

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 6, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                                        16}, ncols, fj, mu, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2042, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2045, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2048, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2051, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2054, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2057, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2060, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2063, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2066, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2069, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2072, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2075, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2078, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2081, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2084, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2087, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2090, 3, 7, 23,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2099, 3, 8, 26,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2108, 3, 9, 29,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2117, 3, 10, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2126, 3, 11, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2135, 3, 12, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2144, 3, 13, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2153, 3, 14, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2162, 3, 15, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2171, 3, 16, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2180, 3, 17, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2189, 3, 18, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2198, 3, 19, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2207, 3, 20, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2216, 3, 21, 65,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2225, 3, 23, 68,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2243, 3, 26, 74,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2261, 3, 29, 80,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2279, 3, 32, 86,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2297, 3, 35, 92,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2315, 3, 38, 98,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2333, 3, 41, 104,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2351, 3, 44, 110,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2369, 3, 47, 116,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2387, 3, 50, 122,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2405, 3, 53, 128,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2423, 3, 56, 134,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2441, 3, 59, 140,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2459, 3, 62, 146,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2477, 3, 68, 152,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2507, 3, 74, 162,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2537, 3, 80, 172,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2567, 3, 86, 182,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2597, 3, 92, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2627, 3, 98, 202,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2657, 3, 104, 212,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2687, 3, 110, 222,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2717, 3, 116, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2747, 3, 122, 242,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2777, 3, 128, 252,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2807, 3, 134, 262,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2837, 3, 140, 272,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2867, 3, 152, 282,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2912, 3, 162, 297,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2957, 3, 172, 312,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3002, 3, 182, 327,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3047, 3, 192, 342,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3092, 3, 202, 357,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3137, 3, 212, 372,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3182, 3, 222, 387,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3227, 3, 232, 402,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3272, 3, 242, 417,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3317, 3, 252, 432,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3362, 3, 262, 447,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3407, 3, 282, 462,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3470, 3, 297, 483,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3533, 3, 312, 504,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3596, 3, 327, 525,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3659, 3, 342, 546,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3722, 3, 357, 567,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3785, 3, 372, 588,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3848, 3, 387, 609,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3911, 3, 402, 630,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3974, 3, 417, 651,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4037, 3, 432, 672,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4100, 3, 462, 693,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4184, 3, 483, 721,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4268, 3, 504, 749,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4352, 3, 525, 777,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4436, 3, 546, 805,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4520, 3, 567, 833,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4604, 3, 588, 861,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4688, 3, 609, 889,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4772, 3, 630, 917,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4856, 3, 651, 945,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4940, 3, 693, 973,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5048, 3, 721,
                                                                       1009, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5156, 3, 749,
                                                                       1045, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5264, 3, 777,
                                                                       1081, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5372, 3, 805,
                                                                       1117, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5480, 3, 833,
                                                                       1153, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5588, 3, 861,
                                                                       1189, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5696, 3, 889,
                                                                       1225, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5804, 3, 917,
                                                                       1261, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5912, 3, 973,
                                                                       1297, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6047, 3, 1009,
                                                                       1342, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6182, 3, 1045,
                                                                       1387, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6317, 3, 1081,
                                                                       1432, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6452, 3, 1117,
                                                                       1477, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6587, 3, 1153,
                                                                       1522, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6722, 3, 1189,
                                                                       1567, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6857, 3, 1225,
                                                                       1612, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 6992, 3, 1297,
                                                                       1657, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7157, 3, 1342,
                                                                       1712, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7322, 3, 1387,
                                                                       1767, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7487, 3, 1432,
                                                                       1822, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7652, 3, 1477,
                                                                       1877, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7817, 3, 1522,
                                                                       1932, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7982, 3, 1567,
                                                                       1987, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8147, 3, 7, 8,
                                                                       2048, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8153, 3, 8, 9,
                                                                       2051, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8159, 3, 9, 10,
                                                                       2054, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8165, 3, 10, 11,
                                                                       2057, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8171, 3, 11, 12,
                                                                       2060, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8177, 3, 12, 13,
                                                                       2063, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8183, 3, 13, 14,
                                                                       2066, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8189, 3, 14, 15,
                                                                       2069, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8195, 3, 15, 16,
                                                                       2072, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8201, 3, 16, 17,
                                                                       2075, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8207, 3, 17, 18,
                                                                       2078, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8213, 3, 18, 19,
                                                                       2081, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8219, 3, 19, 20,
                                                                       2084, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8225, 3, 20, 21,
                                                                       2087, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8231, 0, 3, 8147,
                                                                       2048, 8153, 2108, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8249, 0, 3, 8153,
                                                                       2051, 8159, 2117, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8267, 0, 3, 8159,
                                                                       2054, 8165, 2126, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8285, 0, 3, 8165,
                                                                       2057, 8171, 2135, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8303, 0, 3, 8171,
                                                                       2060, 8177, 2144, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8321, 0, 3, 8177,
                                                                       2063, 8183, 2153, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8339, 0, 3, 8183,
                                                                       2066, 8189, 2162, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8357, 0, 3, 8189,
                                                                       2069, 8195, 2171, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8375, 0, 3, 8195,
                                                                       2072, 8201, 2180, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8393, 0, 3, 8201,
                                                                       2075, 8207, 2189, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8411, 0, 3, 8207,
                                                                       2078, 8213, 2198, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8429, 0, 3, 8213,
                                                                       2081, 8219, 2207, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8447, 0, 3, 8219,
                                                                       2084, 8225, 2216, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8465, 0, 3, 8231,
                                                                       2108, 8249, 68, 74, 2261,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8501, 0, 3, 8249,
                                                                       2117, 8267, 74, 80, 2279,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8537, 0, 3, 8267,
                                                                       2126, 8285, 80, 86, 2297,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8573, 0, 3, 8285,
                                                                       2135, 8303, 86, 92, 2315,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8609, 0, 3, 8303,
                                                                       2144, 8321, 92, 98, 2333,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8645, 0, 3, 8321,
                                                                       2153, 8339, 98, 104, 2351,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8681, 0, 3, 8339,
                                                                       2162, 8357, 104, 110,
                                                                       2369, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8717, 0, 3, 8357,
                                                                       2171, 8375, 110, 116,
                                                                       2387, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8753, 0, 3, 8375,
                                                                       2180, 8393, 116, 122,
                                                                       2405, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8789, 0, 3, 8393,
                                                                       2189, 8411, 122, 128,
                                                                       2423, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8825, 0, 3, 8411,
                                                                       2198, 8429, 128, 134,
                                                                       2441, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8861, 0, 3, 8429,
                                                                       2207, 8447, 134, 140,
                                                                       2459, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8897, 0, 3, 8465,
                                                                       2261, 8501, 152, 162,
                                                                       2537, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8957, 0, 3, 8501,
                                                                       2279, 8537, 162, 172,
                                                                       2567, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9017, 0, 3, 8537,
                                                                       2297, 8573, 172, 182,
                                                                       2597, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9077, 0, 3, 8573,
                                                                       2315, 8609, 182, 192,
                                                                       2627, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9137, 0, 3, 8609,
                                                                       2333, 8645, 192, 202,
                                                                       2657, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9197, 0, 3, 8645,
                                                                       2351, 8681, 202, 212,
                                                                       2687, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9257, 0, 3, 8681,
                                                                       2369, 8717, 212, 222,
                                                                       2717, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9317, 0, 3, 8717,
                                                                       2387, 8753, 222, 232,
                                                                       2747, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9377, 0, 3, 8753,
                                                                       2405, 8789, 232, 242,
                                                                       2777, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9437, 0, 3, 8789,
                                                                       2423, 8825, 242, 252,
                                                                       2807, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9497, 0, 3, 8825,
                                                                       2441, 8861, 252, 262,
                                                                       2837, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9557, 0, 3, 8897,
                                                                       2537, 8957, 282, 297,
                                                                       2957, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9647, 0, 3, 8957,
                                                                       2567, 9017, 297, 312,
                                                                       3002, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9737, 0, 3, 9017,
                                                                       2597, 9077, 312, 327,
                                                                       3047, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9827, 0, 3, 9077,
                                                                       2627, 9137, 327, 342,
                                                                       3092, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9917, 0, 3, 9137,
                                                                       2657, 9197, 342, 357,
                                                                       3137, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10007, 0, 3, 9197,
                                                                       2687, 9257, 357, 372,
                                                                       3182, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10097, 0, 3, 9257,
                                                                       2717, 9317, 372, 387,
                                                                       3227, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10187, 0, 3, 9317,
                                                                       2747, 9377, 387, 402,
                                                                       3272, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10277, 0, 3, 9377,
                                                                       2777, 9437, 402, 417,
                                                                       3317, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10367, 0, 3, 9437,
                                                                       2807, 9497, 417, 432,
                                                                       3362, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10457, 0, 3, 9557,
                                                                       2957, 9647, 462, 483,
                                                                       3533, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10583, 0, 3, 9647,
                                                                       3002, 9737, 483, 504,
                                                                       3596, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10709, 0, 3, 9737,
                                                                       3047, 9827, 504, 525,
                                                                       3659, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10835, 0, 3, 9827,
                                                                       3092, 9917, 525, 546,
                                                                       3722, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10961, 0, 3, 9917,
                                                                       3137, 10007, 546, 567,
                                                                       3785, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11087, 0, 3,
                                                                       10007, 3182, 10097, 567,
                                                                       588, 3848, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11213, 0, 3,
                                                                       10097, 3227, 10187, 588,
                                                                       609, 3911, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11339, 0, 3,
                                                                       10187, 3272, 10277, 609,
                                                                       630, 3974, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11465, 0, 3,
                                                                       10277, 3317, 10367, 630,
                                                                       651, 4037, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11591, 0, 3,
                                                                       10457, 3533, 10583, 693,
                                                                       721, 4268, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11759, 0, 3,
                                                                       10583, 3596, 10709, 721,
                                                                       749, 4352, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11927, 0, 3,
                                                                       10709, 3659, 10835, 749,
                                                                       777, 4436, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12095, 0, 3,
                                                                       10835, 3722, 10961, 777,
                                                                       805, 4520, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12263, 0, 3,
                                                                       10961, 3785, 11087, 805,
                                                                       833, 4604, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12431, 0, 3,
                                                                       11087, 3848, 11213, 833,
                                                                       861, 4688, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12599, 0, 3,
                                                                       11213, 3911, 11339, 861,
                                                                       889, 4772, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12767, 0, 3,
                                                                       11339, 3974, 11465, 889,
                                                                       917, 4856, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12935, 0, 3,
                                                                       11591, 4268, 11759, 973,
                                                                       1009, 5156, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 13151, 0, 3,
                                                                       11759, 4352, 11927, 1009,
                                                                       1045, 5264, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 13367, 0, 3,
                                                                       11927, 4436, 12095, 1045,
                                                                       1081, 5372, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 13583, 0, 3,
                                                                       12095, 4520, 12263, 1081,
                                                                       1117, 5480, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 13799, 0, 3,
                                                                       12263, 4604, 12431, 1117,
                                                                       1153, 5588, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 14015, 0, 3,
                                                                       12431, 4688, 12599, 1153,
                                                                       1189, 5696, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 14231, 0, 3,
                                                                       12599, 4772, 12767, 1189,
                                                                       1225, 5804, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 14447, 0, 3,
                                                                       12935, 5156, 13151, 1297,
                                                                       1342, 6182, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 14717, 0, 3,
                                                                       13151, 5264, 13367, 1342,
                                                                       1387, 6317, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 14987, 0, 3,
                                                                       13367, 5372, 13583, 1387,
                                                                       1432, 6452, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 15257, 0, 3,
                                                                       13583, 5480, 13799, 1432,
                                                                       1477, 6587, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 15527, 0, 3,
                                                                       13799, 5588, 14015, 1477,
                                                                       1522, 6722, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 15797, 0, 3,
                                                                       14015, 5696, 14231, 1522,
                                                                       1567, 6857, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 16067, 0, 3,
                                                                       14447, 6182, 14717, 1657,
                                                                       1712, 7322, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 16397, 0, 3,
                                                                       14717, 6317, 14987, 1712,
                                                                       1767, 7487, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 16727, 0, 3,
                                                                       14987, 6452, 15257, 1767,
                                                                       1822, 7652, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 17057, 0, 3,
                                                                       15257, 6587, 15527, 1822,
                                                                       1877, 7817, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 17387, 0, 3,
                                                                       15527, 6722, 15797, 1877,
                                                                       1932, 7982, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17717, 3, 2042,
                                                                       2045, 8147, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17727, 3, 2045,
                                                                       2048, 8153, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17737, 3, 2048,
                                                                       2051, 8159, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17747, 3, 2051,
                                                                       2054, 8165, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17757, 3, 2054,
                                                                       2057, 8171, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17767, 3, 2057,
                                                                       2060, 8177, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17777, 3, 2060,
                                                                       2063, 8183, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17787, 3, 2063,
                                                                       2066, 8189, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17797, 3, 2066,
                                                                       2069, 8195, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17807, 3, 2069,
                                                                       2072, 8201, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17817, 3, 2072,
                                                                       2075, 8207, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17827, 3, 2075,
                                                                       2078, 8213, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17837, 3, 2078,
                                                                       2081, 8219, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17847, 3, 2081,
                                                                       2084, 8225, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17857, 0, 3,
                                                                       17717, 8147, 17727, 2090,
                                                                       2099, 8231, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17887, 0, 3,
                                                                       17727, 8153, 17737, 2099,
                                                                       2108, 8249, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17917, 0, 3,
                                                                       17737, 8159, 17747, 2108,
                                                                       2117, 8267, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17947, 0, 3,
                                                                       17747, 8165, 17757, 2117,
                                                                       2126, 8285, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17977, 0, 3,
                                                                       17757, 8171, 17767, 2126,
                                                                       2135, 8303, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18007, 0, 3,
                                                                       17767, 8177, 17777, 2135,
                                                                       2144, 8321, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18037, 0, 3,
                                                                       17777, 8183, 17787, 2144,
                                                                       2153, 8339, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18067, 0, 3,
                                                                       17787, 8189, 17797, 2153,
                                                                       2162, 8357, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18097, 0, 3,
                                                                       17797, 8195, 17807, 2162,
                                                                       2171, 8375, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18127, 0, 3,
                                                                       17807, 8201, 17817, 2171,
                                                                       2180, 8393, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18157, 0, 3,
                                                                       17817, 8207, 17827, 2180,
                                                                       2189, 8411, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18187, 0, 3,
                                                                       17827, 8213, 17837, 2189,
                                                                       2198, 8429, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18217, 0, 3,
                                                                       17837, 8219, 17847, 2198,
                                                                       2207, 8447, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18247, 0, 3,
                                                                       17857, 8231, 17887, 2225,
                                                                       2243, 8465, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18307, 0, 3,
                                                                       17887, 8249, 17917, 2243,
                                                                       2261, 8501, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18367, 0, 3,
                                                                       17917, 8267, 17947, 2261,
                                                                       2279, 8537, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18427, 0, 3,
                                                                       17947, 8285, 17977, 2279,
                                                                       2297, 8573, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18487, 0, 3,
                                                                       17977, 8303, 18007, 2297,
                                                                       2315, 8609, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18547, 0, 3,
                                                                       18007, 8321, 18037, 2315,
                                                                       2333, 8645, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18607, 0, 3,
                                                                       18037, 8339, 18067, 2333,
                                                                       2351, 8681, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18667, 0, 3,
                                                                       18067, 8357, 18097, 2351,
                                                                       2369, 8717, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18727, 0, 3,
                                                                       18097, 8375, 18127, 2369,
                                                                       2387, 8753, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18787, 0, 3,
                                                                       18127, 8393, 18157, 2387,
                                                                       2405, 8789, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18847, 0, 3,
                                                                       18157, 8411, 18187, 2405,
                                                                       2423, 8825, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18907, 0, 3,
                                                                       18187, 8429, 18217, 2423,
                                                                       2441, 8861, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18967, 0, 3,
                                                                       18247, 8465, 18307, 2477,
                                                                       2507, 8897, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19067, 0, 3,
                                                                       18307, 8501, 18367, 2507,
                                                                       2537, 8957, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19167, 0, 3,
                                                                       18367, 8537, 18427, 2537,
                                                                       2567, 9017, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19267, 0, 3,
                                                                       18427, 8573, 18487, 2567,
                                                                       2597, 9077, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19367, 0, 3,
                                                                       18487, 8609, 18547, 2597,
                                                                       2627, 9137, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19467, 0, 3,
                                                                       18547, 8645, 18607, 2627,
                                                                       2657, 9197, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19567, 0, 3,
                                                                       18607, 8681, 18667, 2657,
                                                                       2687, 9257, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19667, 0, 3,
                                                                       18667, 8717, 18727, 2687,
                                                                       2717, 9317, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19767, 0, 3,
                                                                       18727, 8753, 18787, 2717,
                                                                       2747, 9377, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19867, 0, 3,
                                                                       18787, 8789, 18847, 2747,
                                                                       2777, 9437, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19967, 0, 3,
                                                                       18847, 8825, 18907, 2777,
                                                                       2807, 9497, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20067, 0, 3,
                                                                       18967, 8897, 19067, 2867,
                                                                       2912, 9557, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20217, 0, 3,
                                                                       19067, 8957, 19167, 2912,
                                                                       2957, 9647, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20367, 0, 3,
                                                                       19167, 9017, 19267, 2957,
                                                                       3002, 9737, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20517, 0, 3,
                                                                       19267, 9077, 19367, 3002,
                                                                       3047, 9827, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20667, 0, 3,
                                                                       19367, 9137, 19467, 3047,
                                                                       3092, 9917, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20817, 0, 3,
                                                                       19467, 9197, 19567, 3092,
                                                                       3137, 10007, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20967, 0, 3,
                                                                       19567, 9257, 19667, 3137,
                                                                       3182, 10097, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 21117, 0, 3,
                                                                       19667, 9317, 19767, 3182,
                                                                       3227, 10187, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 21267, 0, 3,
                                                                       19767, 9377, 19867, 3227,
                                                                       3272, 10277, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 21417, 0, 3,
                                                                       19867, 9437, 19967, 3272,
                                                                       3317, 10367, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 21567, 0, 3,
                                                                       20067, 9557, 20217, 3407,
                                                                       3470, 10457, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 21777, 0, 3,
                                                                       20217, 9647, 20367, 3470,
                                                                       3533, 10583, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 21987, 0, 3,
                                                                       20367, 9737, 20517, 3533,
                                                                       3596, 10709, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 22197, 0, 3,
                                                                       20517, 9827, 20667, 3596,
                                                                       3659, 10835, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 22407, 0, 3,
                                                                       20667, 9917, 20817, 3659,
                                                                       3722, 10961, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 22617, 0, 3,
                                                                       20817, 10007, 20967, 3722,
                                                                       3785, 11087, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 22827, 0, 3,
                                                                       20967, 10097, 21117, 3785,
                                                                       3848, 11213, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 23037, 0, 3,
                                                                       21117, 10187, 21267, 3848,
                                                                       3911, 11339, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 23247, 0, 3,
                                                                       21267, 10277, 21417, 3911,
                                                                       3974, 11465, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 23457, 0, 3,
                                                                       21567, 10457, 21777, 4100,
                                                                       4184, 11591, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 23737, 0, 3,
                                                                       21777, 10583, 21987, 4184,
                                                                       4268, 11759, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 24017, 0, 3,
                                                                       21987, 10709, 22197, 4268,
                                                                       4352, 11927, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 24297, 0, 3,
                                                                       22197, 10835, 22407, 4352,
                                                                       4436, 12095, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 24577, 0, 3,
                                                                       22407, 10961, 22617, 4436,
                                                                       4520, 12263, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 24857, 0, 3,
                                                                       22617, 11087, 22827, 4520,
                                                                       4604, 12431, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 25137, 0, 3,
                                                                       22827, 11213, 23037, 4604,
                                                                       4688, 12599, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 25417, 0, 3,
                                                                       23037, 11339, 23247, 4688,
                                                                       4772, 12767, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 25697, 0, 3,
                                                                       23457, 11591, 23737, 4940,
                                                                       5048, 12935, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 26057, 0, 3,
                                                                       23737, 11759, 24017, 5048,
                                                                       5156, 13151, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 26417, 0, 3,
                                                                       24017, 11927, 24297, 5156,
                                                                       5264, 13367, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 26777, 0, 3,
                                                                       24297, 12095, 24577, 5264,
                                                                       5372, 13583, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 27137, 0, 3,
                                                                       24577, 12263, 24857, 5372,
                                                                       5480, 13799, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 27497, 0, 3,
                                                                       24857, 12431, 25137, 5480,
                                                                       5588, 14015, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 27857, 0, 3,
                                                                       25137, 12599, 25417, 5588,
                                                                       5696, 14231, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 28217, 0, 3,
                                                                       25697, 12935, 26057, 5912,
                                                                       6047, 14447, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 28667, 0, 3,
                                                                       26057, 13151, 26417, 6047,
                                                                       6182, 14717, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 29117, 0, 3,
                                                                       26417, 13367, 26777, 6182,
                                                                       6317, 14987, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 29567, 0, 3,
                                                                       26777, 13583, 27137, 6317,
                                                                       6452, 15257, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 30017, 0, 3,
                                                                       27137, 13799, 27497, 6452,
                                                                       6587, 15527, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 30467, 0, 3,
                                                                       27497, 14015, 27857, 6587,
                                                                       6722, 15797, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 30917, 0, 3,
                                                                       28217, 14447, 28667, 6992,
                                                                       7157, 16067, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 31467, 0, 3,
                                                                       28667, 14717, 29117, 7157,
                                                                       7322, 16397, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 32017, 0, 3,
                                                                       29117, 14987, 29567, 7322,
                                                                       7487, 16727, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 32567, 0, 3,
                                                                       29567, 15257, 30017, 7487,
                                                                       7652, 17057, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 33117, 0, 3,
                                                                       30017, 15527, 30467, 7652,
                                                                       7817, 17387, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33667, 3, 8147,
                                                                       8153, 17737, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33682, 3, 8153,
                                                                       8159, 17747, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33697, 3, 8159,
                                                                       8165, 17757, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33712, 3, 8165,
                                                                       8171, 17767, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33727, 3, 8171,
                                                                       8177, 17777, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33742, 3, 8177,
                                                                       8183, 17787, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33757, 3, 8183,
                                                                       8189, 17797, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33772, 3, 8189,
                                                                       8195, 17807, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33787, 3, 8195,
                                                                       8201, 17817, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33802, 3, 8201,
                                                                       8207, 17827, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33817, 3, 8207,
                                                                       8213, 17837, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33832, 3, 8213,
                                                                       8219, 17847, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 33847, 0, 3,
                                                                       33667, 17737, 33682, 8231,
                                                                       8249, 17917, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 33892, 0, 3,
                                                                       33682, 17747, 33697, 8249,
                                                                       8267, 17947, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 33937, 0, 3,
                                                                       33697, 17757, 33712, 8267,
                                                                       8285, 17977, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 33982, 0, 3,
                                                                       33712, 17767, 33727, 8285,
                                                                       8303, 18007, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34027, 0, 3,
                                                                       33727, 17777, 33742, 8303,
                                                                       8321, 18037, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34072, 0, 3,
                                                                       33742, 17787, 33757, 8321,
                                                                       8339, 18067, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34117, 0, 3,
                                                                       33757, 17797, 33772, 8339,
                                                                       8357, 18097, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34162, 0, 3,
                                                                       33772, 17807, 33787, 8357,
                                                                       8375, 18127, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34207, 0, 3,
                                                                       33787, 17817, 33802, 8375,
                                                                       8393, 18157, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34252, 0, 3,
                                                                       33802, 17827, 33817, 8393,
                                                                       8411, 18187, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34297, 0, 3,
                                                                       33817, 17837, 33832, 8411,
                                                                       8429, 18217, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34342, 0, 3,
                                                                       33847, 17917, 33892, 8465,
                                                                       8501, 18367, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34432, 0, 3,
                                                                       33892, 17947, 33937, 8501,
                                                                       8537, 18427, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34522, 0, 3,
                                                                       33937, 17977, 33982, 8537,
                                                                       8573, 18487, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34612, 0, 3,
                                                                       33982, 18007, 34027, 8573,
                                                                       8609, 18547, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34702, 0, 3,
                                                                       34027, 18037, 34072, 8609,
                                                                       8645, 18607, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34792, 0, 3,
                                                                       34072, 18067, 34117, 8645,
                                                                       8681, 18667, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34882, 0, 3,
                                                                       34117, 18097, 34162, 8681,
                                                                       8717, 18727, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34972, 0, 3,
                                                                       34162, 18127, 34207, 8717,
                                                                       8753, 18787, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 35062, 0, 3,
                                                                       34207, 18157, 34252, 8753,
                                                                       8789, 18847, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 35152, 0, 3,
                                                                       34252, 18187, 34297, 8789,
                                                                       8825, 18907, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 35242, 0, 3,
                                                                       34342, 18367, 34432, 8897,
                                                                       8957, 19167, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 35392, 0, 3,
                                                                       34432, 18427, 34522, 8957,
                                                                       9017, 19267, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 35542, 0, 3,
                                                                       34522, 18487, 34612, 9017,
                                                                       9077, 19367, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 35692, 0, 3,
                                                                       34612, 18547, 34702, 9077,
                                                                       9137, 19467, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 35842, 0, 3,
                                                                       34702, 18607, 34792, 9137,
                                                                       9197, 19567, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 35992, 0, 3,
                                                                       34792, 18667, 34882, 9197,
                                                                       9257, 19667, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 36142, 0, 3,
                                                                       34882, 18727, 34972, 9257,
                                                                       9317, 19767, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 36292, 0, 3,
                                                                       34972, 18787, 35062, 9317,
                                                                       9377, 19867, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 36442, 0, 3,
                                                                       35062, 18847, 35152, 9377,
                                                                       9437, 19967, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 36592, 0, 3,
                                                                       35242, 19167, 35392, 9557,
                                                                       9647, 20367, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 36817, 0, 3,
                                                                       35392, 19267, 35542, 9647,
                                                                       9737, 20517, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 37042, 0, 3,
                                                                       35542, 19367, 35692, 9737,
                                                                       9827, 20667, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 37267, 0, 3,
                                                                       35692, 19467, 35842, 9827,
                                                                       9917, 20817, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 37492, 0, 3,
                                                                       35842, 19567, 35992, 9917,
                                                                       10007, 20967, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 37717, 0, 3,
                                                                       35992, 19667, 36142,
                                                                       10007, 10097, 21117,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 37942, 0, 3,
                                                                       36142, 19767, 36292,
                                                                       10097, 10187, 21267,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 38167, 0, 3,
                                                                       36292, 19867, 36442,
                                                                       10187, 10277, 21417,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 38392, 0, 3,
                                                                       36592, 20367, 36817,
                                                                       10457, 10583, 21987,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 38707, 0, 3,
                                                                       36817, 20517, 37042,
                                                                       10583, 10709, 22197,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 39022, 0, 3,
                                                                       37042, 20667, 37267,
                                                                       10709, 10835, 22407,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 39337, 0, 3,
                                                                       37267, 20817, 37492,
                                                                       10835, 10961, 22617,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 39652, 0, 3,
                                                                       37492, 20967, 37717,
                                                                       10961, 11087, 22827,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 39967, 0, 3,
                                                                       37717, 21117, 37942,
                                                                       11087, 11213, 23037,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 40282, 0, 3,
                                                                       37942, 21267, 38167,
                                                                       11213, 11339, 23247,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 40597, 0, 3,
                                                                       38392, 21987, 38707,
                                                                       11591, 11759, 24017,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 41017, 0, 3,
                                                                       38707, 22197, 39022,
                                                                       11759, 11927, 24297,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 41437, 0, 3,
                                                                       39022, 22407, 39337,
                                                                       11927, 12095, 24577,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 41857, 0, 3,
                                                                       39337, 22617, 39652,
                                                                       12095, 12263, 24857,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 42277, 0, 3,
                                                                       39652, 22827, 39967,
                                                                       12263, 12431, 25137,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 42697, 0, 3,
                                                                       39967, 23037, 40282,
                                                                       12431, 12599, 25417,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 43117, 0, 3,
                                                                       40597, 24017, 41017,
                                                                       12935, 13151, 26417,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 43657, 0, 3,
                                                                       41017, 24297, 41437,
                                                                       13151, 13367, 26777,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 44197, 0, 3,
                                                                       41437, 24577, 41857,
                                                                       13367, 13583, 27137,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 44737, 0, 3,
                                                                       41857, 24857, 42277,
                                                                       13583, 13799, 27497,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 45277, 0, 3,
                                                                       42277, 25137, 42697,
                                                                       13799, 14015, 27857,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 45817, 0, 3,
                                                                       43117, 26417, 43657,
                                                                       14447, 14717, 29117,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 46492, 0, 3,
                                                                       43657, 26777, 44197,
                                                                       14717, 14987, 29567,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 47167, 0, 3,
                                                                       44197, 27137, 44737,
                                                                       14987, 15257, 30017,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 47842, 0, 3,
                                                                       44737, 27497, 45277,
                                                                       15257, 15527, 30467,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 48517, 0, 3,
                                                                       45817, 29117, 46492,
                                                                       16067, 16397, 32017,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 49342, 0, 3,
                                                                       46492, 29567, 47167,
                                                                       16397, 16727, 32567,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 50167, 0, 3,
                                                                       47167, 30017, 47842,
                                                                       16727, 17057, 33117,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50992, 3, 17717,
                                                                       17727, 33667, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51013, 3, 17727,
                                                                       17737, 33682, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51034, 3, 17737,
                                                                       17747, 33697, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51055, 3, 17747,
                                                                       17757, 33712, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51076, 3, 17757,
                                                                       17767, 33727, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51097, 3, 17767,
                                                                       17777, 33742, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51118, 3, 17777,
                                                                       17787, 33757, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51139, 3, 17787,
                                                                       17797, 33772, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51160, 3, 17797,
                                                                       17807, 33787, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51181, 3, 17807,
                                                                       17817, 33802, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51202, 3, 17817,
                                                                       17827, 33817, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51223, 3, 17827,
                                                                       17837, 33832, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51244, 0, 3,
                                                                       50992, 33667, 51013,
                                                                       17857, 17887, 33847,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51307, 0, 3,
                                                                       51013, 33682, 51034,
                                                                       17887, 17917, 33892,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51370, 0, 3,
                                                                       51034, 33697, 51055,
                                                                       17917, 17947, 33937,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51433, 0, 3,
                                                                       51055, 33712, 51076,
                                                                       17947, 17977, 33982,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51496, 0, 3,
                                                                       51076, 33727, 51097,
                                                                       17977, 18007, 34027,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51559, 0, 3,
                                                                       51097, 33742, 51118,
                                                                       18007, 18037, 34072,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51622, 0, 3,
                                                                       51118, 33757, 51139,
                                                                       18037, 18067, 34117,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51685, 0, 3,
                                                                       51139, 33772, 51160,
                                                                       18067, 18097, 34162,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51748, 0, 3,
                                                                       51160, 33787, 51181,
                                                                       18097, 18127, 34207,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51811, 0, 3,
                                                                       51181, 33802, 51202,
                                                                       18127, 18157, 34252,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51874, 0, 3,
                                                                       51202, 33817, 51223,
                                                                       18157, 18187, 34297,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 51937, 0, 3,
                                                                       51244, 33847, 51307,
                                                                       18247, 18307, 34342,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52063, 0, 3,
                                                                       51307, 33892, 51370,
                                                                       18307, 18367, 34432,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52189, 0, 3,
                                                                       51370, 33937, 51433,
                                                                       18367, 18427, 34522,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52315, 0, 3,
                                                                       51433, 33982, 51496,
                                                                       18427, 18487, 34612,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52441, 0, 3,
                                                                       51496, 34027, 51559,
                                                                       18487, 18547, 34702,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52567, 0, 3,
                                                                       51559, 34072, 51622,
                                                                       18547, 18607, 34792,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52693, 0, 3,
                                                                       51622, 34117, 51685,
                                                                       18607, 18667, 34882,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52819, 0, 3,
                                                                       51685, 34162, 51748,
                                                                       18667, 18727, 34972,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52945, 0, 3,
                                                                       51748, 34207, 51811,
                                                                       18727, 18787, 35062,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 53071, 0, 3,
                                                                       51811, 34252, 51874,
                                                                       18787, 18847, 35152,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 53197, 0, 3,
                                                                       51937, 34342, 52063,
                                                                       18967, 19067, 35242,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 53407, 0, 3,
                                                                       52063, 34432, 52189,
                                                                       19067, 19167, 35392,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 53617, 0, 3,
                                                                       52189, 34522, 52315,
                                                                       19167, 19267, 35542,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 53827, 0, 3,
                                                                       52315, 34612, 52441,
                                                                       19267, 19367, 35692,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54037, 0, 3,
                                                                       52441, 34702, 52567,
                                                                       19367, 19467, 35842,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54247, 0, 3,
                                                                       52567, 34792, 52693,
                                                                       19467, 19567, 35992,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54457, 0, 3,
                                                                       52693, 34882, 52819,
                                                                       19567, 19667, 36142,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54667, 0, 3,
                                                                       52819, 34972, 52945,
                                                                       19667, 19767, 36292,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54877, 0, 3,
                                                                       52945, 35062, 53071,
                                                                       19767, 19867, 36442,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 55087, 0, 3,
                                                                       53197, 35242, 53407,
                                                                       20067, 20217, 36592,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 55402, 0, 3,
                                                                       53407, 35392, 53617,
                                                                       20217, 20367, 36817,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 55717, 0, 3,
                                                                       53617, 35542, 53827,
                                                                       20367, 20517, 37042,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 56032, 0, 3,
                                                                       53827, 35692, 54037,
                                                                       20517, 20667, 37267,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 56347, 0, 3,
                                                                       54037, 35842, 54247,
                                                                       20667, 20817, 37492,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 56662, 0, 3,
                                                                       54247, 35992, 54457,
                                                                       20817, 20967, 37717,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 56977, 0, 3,
                                                                       54457, 36142, 54667,
                                                                       20967, 21117, 37942,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 57292, 0, 3,
                                                                       54667, 36292, 54877,
                                                                       21117, 21267, 38167,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 57607, 0, 3,
                                                                       55087, 36592, 55402,
                                                                       21567, 21777, 38392,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 58048, 0, 3,
                                                                       55402, 36817, 55717,
                                                                       21777, 21987, 38707,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 58489, 0, 3,
                                                                       55717, 37042, 56032,
                                                                       21987, 22197, 39022,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 58930, 0, 3,
                                                                       56032, 37267, 56347,
                                                                       22197, 22407, 39337,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 59371, 0, 3,
                                                                       56347, 37492, 56662,
                                                                       22407, 22617, 39652,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 59812, 0, 3,
                                                                       56662, 37717, 56977,
                                                                       22617, 22827, 39967,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 60253, 0, 3,
                                                                       56977, 37942, 57292,
                                                                       22827, 23037, 40282,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 60694, 0, 3,
                                                                       57607, 38392, 58048,
                                                                       23457, 23737, 40597,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 61282, 0, 3,
                                                                       58048, 38707, 58489,
                                                                       23737, 24017, 41017,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 61870, 0, 3,
                                                                       58489, 39022, 58930,
                                                                       24017, 24297, 41437,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 62458, 0, 3,
                                                                       58930, 39337, 59371,
                                                                       24297, 24577, 41857,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 63046, 0, 3,
                                                                       59371, 39652, 59812,
                                                                       24577, 24857, 42277,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 63634, 0, 3,
                                                                       59812, 39967, 60253,
                                                                       24857, 25137, 42697,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 64222, 0, 3,
                                                                       60694, 40597, 61282,
                                                                       25697, 26057, 43117,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 64978, 0, 3,
                                                                       61282, 41017, 61870,
                                                                       26057, 26417, 43657,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 65734, 0, 3,
                                                                       61870, 41437, 62458,
                                                                       26417, 26777, 44197,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 66490, 0, 3,
                                                                       62458, 41857, 63046,
                                                                       26777, 27137, 44737,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 67246, 0, 3,
                                                                       63046, 42277, 63634,
                                                                       27137, 27497, 45277,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 68002, 0, 3,
                                                                       64222, 43117, 64978,
                                                                       28217, 28667, 45817,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 68947, 0, 3,
                                                                       64978, 43657, 65734,
                                                                       28667, 29117, 46492,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 69892, 0, 3,
                                                                       65734, 44197, 66490,
                                                                       29117, 29567, 47167,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 70837, 0, 3,
                                                                       66490, 44737, 67246,
                                                                       29567, 30017, 47842,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 71782, 0, 3,
                                                                       68002, 45817, 68947,
                                                                       30917, 31467, 48517,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 72937, 0, 3,
                                                                       68947, 46492, 69892,
                                                                       31467, 32017, 49342,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 74092, 0, 3,
                                                                       69892, 47167, 70837,
                                                                       32017, 32567, 50167,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75247, 3, 33667,
                                                                       33682, 51034, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75275, 3, 33682,
                                                                       33697, 51055, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75303, 3, 33697,
                                                                       33712, 51076, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75331, 3, 33712,
                                                                       33727, 51097, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75359, 3, 33727,
                                                                       33742, 51118, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75387, 3, 33742,
                                                                       33757, 51139, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75415, 3, 33757,
                                                                       33772, 51160, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75443, 3, 33772,
                                                                       33787, 51181, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75471, 3, 33787,
                                                                       33802, 51202, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75499, 3, 33802,
                                                                       33817, 51223, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 75527, 0, 3,
                                                                       75247, 51034, 75275,
                                                                       33847, 33892, 51370,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 75611, 0, 3,
                                                                       75275, 51055, 75303,
                                                                       33892, 33937, 51433,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 75695, 0, 3,
                                                                       75303, 51076, 75331,
                                                                       33937, 33982, 51496,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 75779, 0, 3,
                                                                       75331, 51097, 75359,
                                                                       33982, 34027, 51559,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 75863, 0, 3,
                                                                       75359, 51118, 75387,
                                                                       34027, 34072, 51622,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 75947, 0, 3,
                                                                       75387, 51139, 75415,
                                                                       34072, 34117, 51685,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 76031, 0, 3,
                                                                       75415, 51160, 75443,
                                                                       34117, 34162, 51748,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 76115, 0, 3,
                                                                       75443, 51181, 75471,
                                                                       34162, 34207, 51811,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 76199, 0, 3,
                                                                       75471, 51202, 75499,
                                                                       34207, 34252, 51874,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 76283, 0, 3,
                                                                       75527, 51370, 75611,
                                                                       34342, 34432, 52189,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 76451, 0, 3,
                                                                       75611, 51433, 75695,
                                                                       34432, 34522, 52315,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 76619, 0, 3,
                                                                       75695, 51496, 75779,
                                                                       34522, 34612, 52441,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 76787, 0, 3,
                                                                       75779, 51559, 75863,
                                                                       34612, 34702, 52567,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 76955, 0, 3,
                                                                       75863, 51622, 75947,
                                                                       34702, 34792, 52693,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 77123, 0, 3,
                                                                       75947, 51685, 76031,
                                                                       34792, 34882, 52819,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 77291, 0, 3,
                                                                       76031, 51748, 76115,
                                                                       34882, 34972, 52945,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 77459, 0, 3,
                                                                       76115, 51811, 76199,
                                                                       34972, 35062, 53071,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 77627, 0, 3,
                                                                       76283, 52189, 76451,
                                                                       35242, 35392, 53617,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 77907, 0, 3,
                                                                       76451, 52315, 76619,
                                                                       35392, 35542, 53827,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 78187, 0, 3,
                                                                       76619, 52441, 76787,
                                                                       35542, 35692, 54037,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 78467, 0, 3,
                                                                       76787, 52567, 76955,
                                                                       35692, 35842, 54247,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 78747, 0, 3,
                                                                       76955, 52693, 77123,
                                                                       35842, 35992, 54457,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 79027, 0, 3,
                                                                       77123, 52819, 77291,
                                                                       35992, 36142, 54667,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 79307, 0, 3,
                                                                       77291, 52945, 77459,
                                                                       36142, 36292, 54877,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 79587, 0, 3,
                                                                       77627, 53617, 77907,
                                                                       36592, 36817, 55717,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 80007, 0, 3,
                                                                       77907, 53827, 78187,
                                                                       36817, 37042, 56032,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 80427, 0, 3,
                                                                       78187, 54037, 78467,
                                                                       37042, 37267, 56347,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 80847, 0, 3,
                                                                       78467, 54247, 78747,
                                                                       37267, 37492, 56662,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 81267, 0, 3,
                                                                       78747, 54457, 79027,
                                                                       37492, 37717, 56977,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 81687, 0, 3,
                                                                       79027, 54667, 79307,
                                                                       37717, 37942, 57292,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 82107, 0, 3,
                                                                       79587, 55717, 80007,
                                                                       38392, 38707, 58489,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 82695, 0, 3,
                                                                       80007, 56032, 80427,
                                                                       38707, 39022, 58930,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 83283, 0, 3,
                                                                       80427, 56347, 80847,
                                                                       39022, 39337, 59371,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 83871, 0, 3,
                                                                       80847, 56662, 81267,
                                                                       39337, 39652, 59812,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 84459, 0, 3,
                                                                       81267, 56977, 81687,
                                                                       39652, 39967, 60253,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 85047, 0, 3,
                                                                       82107, 58489, 82695,
                                                                       40597, 41017, 61870,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 85831, 0, 3,
                                                                       82695, 58930, 83283,
                                                                       41017, 41437, 62458,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 86615, 0, 3,
                                                                       83283, 59371, 83871,
                                                                       41437, 41857, 63046,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 87399, 0, 3,
                                                                       83871, 59812, 84459,
                                                                       41857, 42277, 63634,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 88183, 0, 3,
                                                                       85047, 61870, 85831,
                                                                       43117, 43657, 65734,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 89191, 0, 3,
                                                                       85831, 62458, 86615,
                                                                       43657, 44197, 66490,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 90199, 0, 3,
                                                                       86615, 63046, 87399,
                                                                       44197, 44737, 67246,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 91207, 0, 3,
                                                                       88183, 65734, 89191,
                                                                       45817, 46492, 69892,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 92467, 0, 3,
                                                                       89191, 66490, 90199,
                                                                       46492, 47167, 70837,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 93727, 0, 3,
                                                                       91207, 69892, 92467,
                                                                       48517, 49342, 74092,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95267, 3, 50992,
                                                                       51013, 75247, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95303, 3, 51013,
                                                                       51034, 75275, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95339, 3, 51034,
                                                                       51055, 75303, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95375, 3, 51055,
                                                                       51076, 75331, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95411, 3, 51076,
                                                                       51097, 75359, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95447, 3, 51097,
                                                                       51118, 75387, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95483, 3, 51118,
                                                                       51139, 75415, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95519, 3, 51139,
                                                                       51160, 75443, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95555, 3, 51160,
                                                                       51181, 75471, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95591, 3, 51181,
                                                                       51202, 75499, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 95627, 0, 3,
                                                                       95267, 75247, 95303,
                                                                       51244, 51307, 75527,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 95735, 0, 3,
                                                                       95303, 75275, 95339,
                                                                       51307, 51370, 75611,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 95843, 0, 3,
                                                                       95339, 75303, 95375,
                                                                       51370, 51433, 75695,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 95951, 0, 3,
                                                                       95375, 75331, 95411,
                                                                       51433, 51496, 75779,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 96059, 0, 3,
                                                                       95411, 75359, 95447,
                                                                       51496, 51559, 75863,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 96167, 0, 3,
                                                                       95447, 75387, 95483,
                                                                       51559, 51622, 75947,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 96275, 0, 3,
                                                                       95483, 75415, 95519,
                                                                       51622, 51685, 76031,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 96383, 0, 3,
                                                                       95519, 75443, 95555,
                                                                       51685, 51748, 76115,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 96491, 0, 3,
                                                                       95555, 75471, 95591,
                                                                       51748, 51811, 76199,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 96599, 0, 3,
                                                                       95627, 75527, 95735,
                                                                       51937, 52063, 76283,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 96815, 0, 3,
                                                                       95735, 75611, 95843,
                                                                       52063, 52189, 76451,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 97031, 0, 3,
                                                                       95843, 75695, 95951,
                                                                       52189, 52315, 76619,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 97247, 0, 3,
                                                                       95951, 75779, 96059,
                                                                       52315, 52441, 76787,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 97463, 0, 3,
                                                                       96059, 75863, 96167,
                                                                       52441, 52567, 76955,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 97679, 0, 3,
                                                                       96167, 75947, 96275,
                                                                       52567, 52693, 77123,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 97895, 0, 3,
                                                                       96275, 76031, 96383,
                                                                       52693, 52819, 77291,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 98111, 0, 3,
                                                                       96383, 76115, 96491,
                                                                       52819, 52945, 77459,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 98327, 0, 3,
                                                                       96599, 76283, 96815,
                                                                       53197, 53407, 77627,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 98687, 0, 3,
                                                                       96815, 76451, 97031,
                                                                       53407, 53617, 77907,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 99047, 0, 3,
                                                                       97031, 76619, 97247,
                                                                       53617, 53827, 78187,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 99407, 0, 3,
                                                                       97247, 76787, 97463,
                                                                       53827, 54037, 78467,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 99767, 0, 3,
                                                                       97463, 76955, 97679,
                                                                       54037, 54247, 78747,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 100127, 0, 3,
                                                                       97679, 77123, 97895,
                                                                       54247, 54457, 79027,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 100487, 0, 3,
                                                                       97895, 77291, 98111,
                                                                       54457, 54667, 79307,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 100847, 0, 3,
                                                                       98327, 77627, 98687,
                                                                       55087, 55402, 79587,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 101387, 0, 3,
                                                                       98687, 77907, 99047,
                                                                       55402, 55717, 80007,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 101927, 0, 3,
                                                                       99047, 78187, 99407,
                                                                       55717, 56032, 80427,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 102467, 0, 3,
                                                                       99407, 78467, 99767,
                                                                       56032, 56347, 80847,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 103007, 0, 3,
                                                                       99767, 78747, 100127,
                                                                       56347, 56662, 81267,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 103547, 0, 3,
                                                                       100127, 79027, 100487,
                                                                       56662, 56977, 81687,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 104087, 0, 3,
                                                                       100847, 79587, 101387,
                                                                       57607, 58048, 82107,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 104843, 0, 3,
                                                                       101387, 80007, 101927,
                                                                       58048, 58489, 82695,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 105599, 0, 3,
                                                                       101927, 80427, 102467,
                                                                       58489, 58930, 83283,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 106355, 0, 3,
                                                                       102467, 80847, 103007,
                                                                       58930, 59371, 83871,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 107111, 0, 3,
                                                                       103007, 81267, 103547,
                                                                       59371, 59812, 84459,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 107867, 0, 3,
                                                                       104087, 82107, 104843,
                                                                       60694, 61282, 85047,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 108875, 0, 3,
                                                                       104843, 82695, 105599,
                                                                       61282, 61870, 85831,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 109883, 0, 3,
                                                                       105599, 83283, 106355,
                                                                       61870, 62458, 86615,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 110891, 0, 3,
                                                                       106355, 83871, 107111,
                                                                       62458, 63046, 87399,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 111899, 0, 3,
                                                                       107867, 85047, 108875,
                                                                       64222, 64978, 88183,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 113195, 0, 3,
                                                                       108875, 85831, 109883,
                                                                       64978, 65734, 89191,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 114491, 0, 3,
                                                                       109883, 86615, 110891,
                                                                       65734, 66490, 90199,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 115787, 0, 3,
                                                                       111899, 88183, 113195,
                                                                       68002, 68947, 91207,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 117407, 0, 3,
                                                                       113195, 89191, 114491,
                                                                       68947, 69892, 92467,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 119027, 0, 3,
                                                                       115787, 91207, 117407,
                                                                       71782, 72937, 93727,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 121007, 107867, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 122435, 111899, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 124271, 115787, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 126566, 119027, 1980, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 122015, 121007, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 123731, 122435, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 125891, 124271, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 128546, 126566, 55, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 129371, 122015, 123731, 15,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 130631, 123731, 125891, 15,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 132251, 125891, 128546, 15,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 134276, 129371, 130631, 15,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 136796, 130631, 132251, 15,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 140036, 134276, 136796, 15,
                                             nmax);

        simdtrf::transform_f_inner(buffer, 144236, 140036, 28, 15, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 144236, 105, nmax);
    }

    for (size_t m = 0; m < 1365; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
