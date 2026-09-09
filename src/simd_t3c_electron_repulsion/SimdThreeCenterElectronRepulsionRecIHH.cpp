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


#include "SimdThreeCenterElectronRepulsionRecIHH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
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
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_ihh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_ihh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 154595, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1573 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 154595, 93919, 8998, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2828, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2831, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2834, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2837, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2840, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2843, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2846, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2849, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2852, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2855, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2858, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2861, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2864, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2867, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2870, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2873, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2876, 3, 7, 23,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2885, 3, 8, 26,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2894, 3, 9, 29,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2903, 3, 10, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2912, 3, 11, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2921, 3, 12, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2930, 3, 13, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2939, 3, 14, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2948, 3, 15, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2957, 3, 16, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2966, 3, 17, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2975, 3, 18, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2984, 3, 19, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2993, 3, 20, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3002, 3, 21, 65,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3011, 3, 23, 68,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3029, 3, 26, 74,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3047, 3, 29, 80,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3065, 3, 32, 86,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3083, 3, 35, 92,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3101, 3, 38, 98,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3119, 3, 41, 104,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3137, 3, 44, 110,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3155, 3, 47, 116,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3173, 3, 50, 122,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3191, 3, 53, 128,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3209, 3, 56, 134,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3227, 3, 59, 140,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3245, 3, 62, 146,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3263, 3, 68, 152,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3293, 3, 74, 162,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3323, 3, 80, 172,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3353, 3, 86, 182,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3383, 3, 92, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3413, 3, 98, 202,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3443, 3, 104, 212,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3473, 3, 110, 222,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3503, 3, 116, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3533, 3, 122, 242,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3563, 3, 128, 252,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3593, 3, 134, 262,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3623, 3, 140, 272,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3653, 3, 152, 282,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3698, 3, 162, 297,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3743, 3, 172, 312,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3788, 3, 182, 327,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3833, 3, 192, 342,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3878, 3, 202, 357,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3923, 3, 212, 372,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3968, 3, 222, 387,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4013, 3, 232, 402,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4058, 3, 242, 417,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4103, 3, 252, 432,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4148, 3, 262, 447,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4193, 3, 282, 462,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4256, 3, 297, 483,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4319, 3, 312, 504,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4382, 3, 327, 525,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4445, 3, 342, 546,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4508, 3, 357, 567,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4571, 3, 372, 588,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4634, 3, 387, 609,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4697, 3, 402, 630,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4760, 3, 417, 651,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4823, 3, 432, 672,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4886, 3, 462, 693,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4970, 3, 483, 721,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5054, 3, 504, 749,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5138, 3, 525, 777,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5222, 3, 546, 805,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5306, 3, 567, 833,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5390, 3, 588, 861,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5474, 3, 609, 889,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5558, 3, 630, 917,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5642, 3, 651, 945,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5726, 3, 693, 973,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5834, 3, 721,
                                                                       1009, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5942, 3, 749,
                                                                       1045, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6050, 3, 777,
                                                                       1081, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6158, 3, 805,
                                                                       1117, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6266, 3, 833,
                                                                       1153, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6374, 3, 861,
                                                                       1189, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6482, 3, 889,
                                                                       1225, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6590, 3, 917,
                                                                       1261, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6698, 3, 973,
                                                                       1297, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6833, 3, 1009,
                                                                       1342, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6968, 3, 1045,
                                                                       1387, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7103, 3, 1081,
                                                                       1432, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7238, 3, 1117,
                                                                       1477, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7373, 3, 1153,
                                                                       1522, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7508, 3, 1189,
                                                                       1567, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7643, 3, 1225,
                                                                       1612, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7778, 3, 1297,
                                                                       1657, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7943, 3, 1342,
                                                                       1712, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8108, 3, 1387,
                                                                       1767, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8273, 3, 1432,
                                                                       1822, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8438, 3, 1477,
                                                                       1877, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8603, 3, 1522,
                                                                       1932, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8768, 3, 1567,
                                                                       1987, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 8933, 3, 1657,
                                                                       2042, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9131, 3, 1712,
                                                                       2108, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9329, 3, 1767,
                                                                       2174, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9527, 3, 1822,
                                                                       2240, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9725, 3, 1877,
                                                                       2306, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 9923, 3, 1932,
                                                                       2372, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 10121, 3, 2042,
                                                                       2438, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 10355, 3, 2108,
                                                                       2516, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 10589, 3, 2174,
                                                                       2594, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 10823, 3, 2240,
                                                                       2672, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 11057, 3, 2306,
                                                                       2750, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11291, 3, 7, 8,
                                                                       2834, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11297, 3, 8, 9,
                                                                       2837, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11303, 3, 9, 10,
                                                                       2840, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11309, 3, 10, 11,
                                                                       2843, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11315, 3, 11, 12,
                                                                       2846, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11321, 3, 12, 13,
                                                                       2849, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11327, 3, 13, 14,
                                                                       2852, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11333, 3, 14, 15,
                                                                       2855, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11339, 3, 15, 16,
                                                                       2858, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11345, 3, 16, 17,
                                                                       2861, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11351, 3, 17, 18,
                                                                       2864, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11357, 3, 18, 19,
                                                                       2867, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11363, 3, 19, 20,
                                                                       2870, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11369, 3, 20, 21,
                                                                       2873, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11375, 0, 3,
                                                                       11291, 2834, 11297, 2894,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11393, 0, 3,
                                                                       11297, 2837, 11303, 2903,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11411, 0, 3,
                                                                       11303, 2840, 11309, 2912,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11429, 0, 3,
                                                                       11309, 2843, 11315, 2921,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11447, 0, 3,
                                                                       11315, 2846, 11321, 2930,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11465, 0, 3,
                                                                       11321, 2849, 11327, 2939,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11483, 0, 3,
                                                                       11327, 2852, 11333, 2948,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11501, 0, 3,
                                                                       11333, 2855, 11339, 2957,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11519, 0, 3,
                                                                       11339, 2858, 11345, 2966,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11537, 0, 3,
                                                                       11345, 2861, 11351, 2975,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11555, 0, 3,
                                                                       11351, 2864, 11357, 2984,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11573, 0, 3,
                                                                       11357, 2867, 11363, 2993,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11591, 0, 3,
                                                                       11363, 2870, 11369, 3002,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11609, 0, 3,
                                                                       11375, 2894, 11393, 68,
                                                                       74, 3047, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11645, 0, 3,
                                                                       11393, 2903, 11411, 74,
                                                                       80, 3065, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11681, 0, 3,
                                                                       11411, 2912, 11429, 80,
                                                                       86, 3083, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11717, 0, 3,
                                                                       11429, 2921, 11447, 86,
                                                                       92, 3101, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11753, 0, 3,
                                                                       11447, 2930, 11465, 92,
                                                                       98, 3119, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11789, 0, 3,
                                                                       11465, 2939, 11483, 98,
                                                                       104, 3137, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11825, 0, 3,
                                                                       11483, 2948, 11501, 104,
                                                                       110, 3155, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11861, 0, 3,
                                                                       11501, 2957, 11519, 110,
                                                                       116, 3173, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11897, 0, 3,
                                                                       11519, 2966, 11537, 116,
                                                                       122, 3191, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11933, 0, 3,
                                                                       11537, 2975, 11555, 122,
                                                                       128, 3209, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11969, 0, 3,
                                                                       11555, 2984, 11573, 128,
                                                                       134, 3227, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12005, 0, 3,
                                                                       11573, 2993, 11591, 134,
                                                                       140, 3245, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12041, 0, 3,
                                                                       11609, 3047, 11645, 152,
                                                                       162, 3323, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12101, 0, 3,
                                                                       11645, 3065, 11681, 162,
                                                                       172, 3353, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12161, 0, 3,
                                                                       11681, 3083, 11717, 172,
                                                                       182, 3383, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12221, 0, 3,
                                                                       11717, 3101, 11753, 182,
                                                                       192, 3413, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12281, 0, 3,
                                                                       11753, 3119, 11789, 192,
                                                                       202, 3443, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12341, 0, 3,
                                                                       11789, 3137, 11825, 202,
                                                                       212, 3473, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12401, 0, 3,
                                                                       11825, 3155, 11861, 212,
                                                                       222, 3503, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12461, 0, 3,
                                                                       11861, 3173, 11897, 222,
                                                                       232, 3533, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12521, 0, 3,
                                                                       11897, 3191, 11933, 232,
                                                                       242, 3563, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12581, 0, 3,
                                                                       11933, 3209, 11969, 242,
                                                                       252, 3593, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12641, 0, 3,
                                                                       11969, 3227, 12005, 252,
                                                                       262, 3623, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12701, 0, 3,
                                                                       12041, 3323, 12101, 282,
                                                                       297, 3743, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12791, 0, 3,
                                                                       12101, 3353, 12161, 297,
                                                                       312, 3788, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12881, 0, 3,
                                                                       12161, 3383, 12221, 312,
                                                                       327, 3833, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12971, 0, 3,
                                                                       12221, 3413, 12281, 327,
                                                                       342, 3878, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13061, 0, 3,
                                                                       12281, 3443, 12341, 342,
                                                                       357, 3923, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13151, 0, 3,
                                                                       12341, 3473, 12401, 357,
                                                                       372, 3968, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13241, 0, 3,
                                                                       12401, 3503, 12461, 372,
                                                                       387, 4013, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13331, 0, 3,
                                                                       12461, 3533, 12521, 387,
                                                                       402, 4058, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13421, 0, 3,
                                                                       12521, 3563, 12581, 402,
                                                                       417, 4103, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13511, 0, 3,
                                                                       12581, 3593, 12641, 417,
                                                                       432, 4148, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13601, 0, 3,
                                                                       12701, 3743, 12791, 462,
                                                                       483, 4319, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13727, 0, 3,
                                                                       12791, 3788, 12881, 483,
                                                                       504, 4382, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13853, 0, 3,
                                                                       12881, 3833, 12971, 504,
                                                                       525, 4445, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13979, 0, 3,
                                                                       12971, 3878, 13061, 525,
                                                                       546, 4508, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14105, 0, 3,
                                                                       13061, 3923, 13151, 546,
                                                                       567, 4571, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14231, 0, 3,
                                                                       13151, 3968, 13241, 567,
                                                                       588, 4634, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14357, 0, 3,
                                                                       13241, 4013, 13331, 588,
                                                                       609, 4697, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14483, 0, 3,
                                                                       13331, 4058, 13421, 609,
                                                                       630, 4760, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14609, 0, 3,
                                                                       13421, 4103, 13511, 630,
                                                                       651, 4823, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14735, 0, 3,
                                                                       13601, 4319, 13727, 693,
                                                                       721, 5054, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14903, 0, 3,
                                                                       13727, 4382, 13853, 721,
                                                                       749, 5138, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15071, 0, 3,
                                                                       13853, 4445, 13979, 749,
                                                                       777, 5222, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15239, 0, 3,
                                                                       13979, 4508, 14105, 777,
                                                                       805, 5306, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15407, 0, 3,
                                                                       14105, 4571, 14231, 805,
                                                                       833, 5390, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15575, 0, 3,
                                                                       14231, 4634, 14357, 833,
                                                                       861, 5474, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15743, 0, 3,
                                                                       14357, 4697, 14483, 861,
                                                                       889, 5558, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15911, 0, 3,
                                                                       14483, 4760, 14609, 889,
                                                                       917, 5642, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16079, 0, 3,
                                                                       14735, 5054, 14903, 973,
                                                                       1009, 5942, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16295, 0, 3,
                                                                       14903, 5138, 15071, 1009,
                                                                       1045, 6050, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16511, 0, 3,
                                                                       15071, 5222, 15239, 1045,
                                                                       1081, 6158, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16727, 0, 3,
                                                                       15239, 5306, 15407, 1081,
                                                                       1117, 6266, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16943, 0, 3,
                                                                       15407, 5390, 15575, 1117,
                                                                       1153, 6374, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17159, 0, 3,
                                                                       15575, 5474, 15743, 1153,
                                                                       1189, 6482, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17375, 0, 3,
                                                                       15743, 5558, 15911, 1189,
                                                                       1225, 6590, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17591, 0, 3,
                                                                       16079, 5942, 16295, 1297,
                                                                       1342, 6968, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17861, 0, 3,
                                                                       16295, 6050, 16511, 1342,
                                                                       1387, 7103, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18131, 0, 3,
                                                                       16511, 6158, 16727, 1387,
                                                                       1432, 7238, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18401, 0, 3,
                                                                       16727, 6266, 16943, 1432,
                                                                       1477, 7373, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18671, 0, 3,
                                                                       16943, 6374, 17159, 1477,
                                                                       1522, 7508, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18941, 0, 3,
                                                                       17159, 6482, 17375, 1522,
                                                                       1567, 7643, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 19211, 0, 3,
                                                                       17591, 6968, 17861, 1657,
                                                                       1712, 8108, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 19541, 0, 3,
                                                                       17861, 7103, 18131, 1712,
                                                                       1767, 8273, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 19871, 0, 3,
                                                                       18131, 7238, 18401, 1767,
                                                                       1822, 8438, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 20201, 0, 3,
                                                                       18401, 7373, 18671, 1822,
                                                                       1877, 8603, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 20531, 0, 3,
                                                                       18671, 7508, 18941, 1877,
                                                                       1932, 8768, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 20861, 0, 3,
                                                                       19211, 8108, 19541, 2042,
                                                                       2108, 9329, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 21257, 0, 3,
                                                                       19541, 8273, 19871, 2108,
                                                                       2174, 9527, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 21653, 0, 3,
                                                                       19871, 8438, 20201, 2174,
                                                                       2240, 9725, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 22049, 0, 3,
                                                                       20201, 8603, 20531, 2240,
                                                                       2306, 9923, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 22445, 0, 3,
                                                                       20861, 9329, 21257, 2438,
                                                                       2516, 10589, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 22913, 0, 3,
                                                                       21257, 9527, 21653, 2516,
                                                                       2594, 10823, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 23381, 0, 3,
                                                                       21653, 9725, 22049, 2594,
                                                                       2672, 11057, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23849, 3, 2828,
                                                                       2831, 11291, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23859, 3, 2831,
                                                                       2834, 11297, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23869, 3, 2834,
                                                                       2837, 11303, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23879, 3, 2837,
                                                                       2840, 11309, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23889, 3, 2840,
                                                                       2843, 11315, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23899, 3, 2843,
                                                                       2846, 11321, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23909, 3, 2846,
                                                                       2849, 11327, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23919, 3, 2849,
                                                                       2852, 11333, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23929, 3, 2852,
                                                                       2855, 11339, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23939, 3, 2855,
                                                                       2858, 11345, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23949, 3, 2858,
                                                                       2861, 11351, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23959, 3, 2861,
                                                                       2864, 11357, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23969, 3, 2864,
                                                                       2867, 11363, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23979, 3, 2867,
                                                                       2870, 11369, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 23989, 0, 3,
                                                                       23849, 11291, 23859, 2876,
                                                                       2885, 11375, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24019, 0, 3,
                                                                       23859, 11297, 23869, 2885,
                                                                       2894, 11393, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24049, 0, 3,
                                                                       23869, 11303, 23879, 2894,
                                                                       2903, 11411, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24079, 0, 3,
                                                                       23879, 11309, 23889, 2903,
                                                                       2912, 11429, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24109, 0, 3,
                                                                       23889, 11315, 23899, 2912,
                                                                       2921, 11447, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24139, 0, 3,
                                                                       23899, 11321, 23909, 2921,
                                                                       2930, 11465, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24169, 0, 3,
                                                                       23909, 11327, 23919, 2930,
                                                                       2939, 11483, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24199, 0, 3,
                                                                       23919, 11333, 23929, 2939,
                                                                       2948, 11501, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24229, 0, 3,
                                                                       23929, 11339, 23939, 2948,
                                                                       2957, 11519, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24259, 0, 3,
                                                                       23939, 11345, 23949, 2957,
                                                                       2966, 11537, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24289, 0, 3,
                                                                       23949, 11351, 23959, 2966,
                                                                       2975, 11555, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24319, 0, 3,
                                                                       23959, 11357, 23969, 2975,
                                                                       2984, 11573, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24349, 0, 3,
                                                                       23969, 11363, 23979, 2984,
                                                                       2993, 11591, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24379, 0, 3,
                                                                       23989, 11375, 24019, 3011,
                                                                       3029, 11609, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24439, 0, 3,
                                                                       24019, 11393, 24049, 3029,
                                                                       3047, 11645, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24499, 0, 3,
                                                                       24049, 11411, 24079, 3047,
                                                                       3065, 11681, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24559, 0, 3,
                                                                       24079, 11429, 24109, 3065,
                                                                       3083, 11717, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24619, 0, 3,
                                                                       24109, 11447, 24139, 3083,
                                                                       3101, 11753, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24679, 0, 3,
                                                                       24139, 11465, 24169, 3101,
                                                                       3119, 11789, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24739, 0, 3,
                                                                       24169, 11483, 24199, 3119,
                                                                       3137, 11825, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24799, 0, 3,
                                                                       24199, 11501, 24229, 3137,
                                                                       3155, 11861, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24859, 0, 3,
                                                                       24229, 11519, 24259, 3155,
                                                                       3173, 11897, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24919, 0, 3,
                                                                       24259, 11537, 24289, 3173,
                                                                       3191, 11933, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24979, 0, 3,
                                                                       24289, 11555, 24319, 3191,
                                                                       3209, 11969, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 25039, 0, 3,
                                                                       24319, 11573, 24349, 3209,
                                                                       3227, 12005, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25099, 0, 3,
                                                                       24379, 11609, 24439, 3263,
                                                                       3293, 12041, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25199, 0, 3,
                                                                       24439, 11645, 24499, 3293,
                                                                       3323, 12101, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25299, 0, 3,
                                                                       24499, 11681, 24559, 3323,
                                                                       3353, 12161, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25399, 0, 3,
                                                                       24559, 11717, 24619, 3353,
                                                                       3383, 12221, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25499, 0, 3,
                                                                       24619, 11753, 24679, 3383,
                                                                       3413, 12281, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25599, 0, 3,
                                                                       24679, 11789, 24739, 3413,
                                                                       3443, 12341, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25699, 0, 3,
                                                                       24739, 11825, 24799, 3443,
                                                                       3473, 12401, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25799, 0, 3,
                                                                       24799, 11861, 24859, 3473,
                                                                       3503, 12461, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25899, 0, 3,
                                                                       24859, 11897, 24919, 3503,
                                                                       3533, 12521, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25999, 0, 3,
                                                                       24919, 11933, 24979, 3533,
                                                                       3563, 12581, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 26099, 0, 3,
                                                                       24979, 11969, 25039, 3563,
                                                                       3593, 12641, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26199, 0, 3,
                                                                       25099, 12041, 25199, 3653,
                                                                       3698, 12701, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26349, 0, 3,
                                                                       25199, 12101, 25299, 3698,
                                                                       3743, 12791, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26499, 0, 3,
                                                                       25299, 12161, 25399, 3743,
                                                                       3788, 12881, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26649, 0, 3,
                                                                       25399, 12221, 25499, 3788,
                                                                       3833, 12971, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26799, 0, 3,
                                                                       25499, 12281, 25599, 3833,
                                                                       3878, 13061, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26949, 0, 3,
                                                                       25599, 12341, 25699, 3878,
                                                                       3923, 13151, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27099, 0, 3,
                                                                       25699, 12401, 25799, 3923,
                                                                       3968, 13241, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27249, 0, 3,
                                                                       25799, 12461, 25899, 3968,
                                                                       4013, 13331, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27399, 0, 3,
                                                                       25899, 12521, 25999, 4013,
                                                                       4058, 13421, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27549, 0, 3,
                                                                       25999, 12581, 26099, 4058,
                                                                       4103, 13511, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 27699, 0, 3,
                                                                       26199, 12701, 26349, 4193,
                                                                       4256, 13601, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 27909, 0, 3,
                                                                       26349, 12791, 26499, 4256,
                                                                       4319, 13727, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 28119, 0, 3,
                                                                       26499, 12881, 26649, 4319,
                                                                       4382, 13853, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 28329, 0, 3,
                                                                       26649, 12971, 26799, 4382,
                                                                       4445, 13979, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 28539, 0, 3,
                                                                       26799, 13061, 26949, 4445,
                                                                       4508, 14105, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 28749, 0, 3,
                                                                       26949, 13151, 27099, 4508,
                                                                       4571, 14231, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 28959, 0, 3,
                                                                       27099, 13241, 27249, 4571,
                                                                       4634, 14357, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 29169, 0, 3,
                                                                       27249, 13331, 27399, 4634,
                                                                       4697, 14483, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 29379, 0, 3,
                                                                       27399, 13421, 27549, 4697,
                                                                       4760, 14609, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 29589, 0, 3,
                                                                       27699, 13601, 27909, 4886,
                                                                       4970, 14735, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 29869, 0, 3,
                                                                       27909, 13727, 28119, 4970,
                                                                       5054, 14903, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 30149, 0, 3,
                                                                       28119, 13853, 28329, 5054,
                                                                       5138, 15071, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 30429, 0, 3,
                                                                       28329, 13979, 28539, 5138,
                                                                       5222, 15239, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 30709, 0, 3,
                                                                       28539, 14105, 28749, 5222,
                                                                       5306, 15407, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 30989, 0, 3,
                                                                       28749, 14231, 28959, 5306,
                                                                       5390, 15575, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 31269, 0, 3,
                                                                       28959, 14357, 29169, 5390,
                                                                       5474, 15743, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 31549, 0, 3,
                                                                       29169, 14483, 29379, 5474,
                                                                       5558, 15911, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 31829, 0, 3,
                                                                       29589, 14735, 29869, 5726,
                                                                       5834, 16079, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 32189, 0, 3,
                                                                       29869, 14903, 30149, 5834,
                                                                       5942, 16295, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 32549, 0, 3,
                                                                       30149, 15071, 30429, 5942,
                                                                       6050, 16511, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 32909, 0, 3,
                                                                       30429, 15239, 30709, 6050,
                                                                       6158, 16727, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 33269, 0, 3,
                                                                       30709, 15407, 30989, 6158,
                                                                       6266, 16943, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 33629, 0, 3,
                                                                       30989, 15575, 31269, 6266,
                                                                       6374, 17159, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 33989, 0, 3,
                                                                       31269, 15743, 31549, 6374,
                                                                       6482, 17375, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 34349, 0, 3,
                                                                       31829, 16079, 32189, 6698,
                                                                       6833, 17591, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 34799, 0, 3,
                                                                       32189, 16295, 32549, 6833,
                                                                       6968, 17861, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 35249, 0, 3,
                                                                       32549, 16511, 32909, 6968,
                                                                       7103, 18131, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 35699, 0, 3,
                                                                       32909, 16727, 33269, 7103,
                                                                       7238, 18401, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 36149, 0, 3,
                                                                       33269, 16943, 33629, 7238,
                                                                       7373, 18671, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 36599, 0, 3,
                                                                       33629, 17159, 33989, 7373,
                                                                       7508, 18941, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 37049, 0, 3,
                                                                       34349, 17591, 34799, 7778,
                                                                       7943, 19211, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 37599, 0, 3,
                                                                       34799, 17861, 35249, 7943,
                                                                       8108, 19541, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 38149, 0, 3,
                                                                       35249, 18131, 35699, 8108,
                                                                       8273, 19871, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 38699, 0, 3,
                                                                       35699, 18401, 36149, 8273,
                                                                       8438, 20201, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 39249, 0, 3,
                                                                       36149, 18671, 36599, 8438,
                                                                       8603, 20531, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 39799, 0, 3,
                                                                       37049, 19211, 37599, 8933,
                                                                       9131, 20861, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 40459, 0, 3,
                                                                       37599, 19541, 38149, 9131,
                                                                       9329, 21257, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 41119, 0, 3,
                                                                       38149, 19871, 38699, 9329,
                                                                       9527, 21653, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 41779, 0, 3,
                                                                       38699, 20201, 39249, 9527,
                                                                       9725, 22049, ncols, gamma,
                                                                       p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 42439, 0, 3,
                                                                       39799, 20861, 40459,
                                                                       10121, 10355, 22445,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 43219, 0, 3,
                                                                       40459, 21257, 41119,
                                                                       10355, 10589, 22913,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 43999, 0, 3,
                                                                       41119, 21653, 41779,
                                                                       10589, 10823, 23381,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44779, 3, 11291,
                                                                       11297, 23869, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44794, 3, 11297,
                                                                       11303, 23879, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44809, 3, 11303,
                                                                       11309, 23889, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44824, 3, 11309,
                                                                       11315, 23899, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44839, 3, 11315,
                                                                       11321, 23909, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44854, 3, 11321,
                                                                       11327, 23919, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44869, 3, 11327,
                                                                       11333, 23929, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44884, 3, 11333,
                                                                       11339, 23939, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44899, 3, 11339,
                                                                       11345, 23949, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44914, 3, 11345,
                                                                       11351, 23959, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44929, 3, 11351,
                                                                       11357, 23969, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44944, 3, 11357,
                                                                       11363, 23979, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 44959, 0, 3,
                                                                       44779, 23869, 44794,
                                                                       11375, 11393, 24049,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45004, 0, 3,
                                                                       44794, 23879, 44809,
                                                                       11393, 11411, 24079,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45049, 0, 3,
                                                                       44809, 23889, 44824,
                                                                       11411, 11429, 24109,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45094, 0, 3,
                                                                       44824, 23899, 44839,
                                                                       11429, 11447, 24139,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45139, 0, 3,
                                                                       44839, 23909, 44854,
                                                                       11447, 11465, 24169,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45184, 0, 3,
                                                                       44854, 23919, 44869,
                                                                       11465, 11483, 24199,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45229, 0, 3,
                                                                       44869, 23929, 44884,
                                                                       11483, 11501, 24229,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45274, 0, 3,
                                                                       44884, 23939, 44899,
                                                                       11501, 11519, 24259,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45319, 0, 3,
                                                                       44899, 23949, 44914,
                                                                       11519, 11537, 24289,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45364, 0, 3,
                                                                       44914, 23959, 44929,
                                                                       11537, 11555, 24319,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 45409, 0, 3,
                                                                       44929, 23969, 44944,
                                                                       11555, 11573, 24349,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 45454, 0, 3,
                                                                       44959, 24049, 45004,
                                                                       11609, 11645, 24499,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 45544, 0, 3,
                                                                       45004, 24079, 45049,
                                                                       11645, 11681, 24559,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 45634, 0, 3,
                                                                       45049, 24109, 45094,
                                                                       11681, 11717, 24619,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 45724, 0, 3,
                                                                       45094, 24139, 45139,
                                                                       11717, 11753, 24679,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 45814, 0, 3,
                                                                       45139, 24169, 45184,
                                                                       11753, 11789, 24739,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 45904, 0, 3,
                                                                       45184, 24199, 45229,
                                                                       11789, 11825, 24799,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 45994, 0, 3,
                                                                       45229, 24229, 45274,
                                                                       11825, 11861, 24859,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46084, 0, 3,
                                                                       45274, 24259, 45319,
                                                                       11861, 11897, 24919,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46174, 0, 3,
                                                                       45319, 24289, 45364,
                                                                       11897, 11933, 24979,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 46264, 0, 3,
                                                                       45364, 24319, 45409,
                                                                       11933, 11969, 25039,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 46354, 0, 3,
                                                                       45454, 24499, 45544,
                                                                       12041, 12101, 25299,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 46504, 0, 3,
                                                                       45544, 24559, 45634,
                                                                       12101, 12161, 25399,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 46654, 0, 3,
                                                                       45634, 24619, 45724,
                                                                       12161, 12221, 25499,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 46804, 0, 3,
                                                                       45724, 24679, 45814,
                                                                       12221, 12281, 25599,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 46954, 0, 3,
                                                                       45814, 24739, 45904,
                                                                       12281, 12341, 25699,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47104, 0, 3,
                                                                       45904, 24799, 45994,
                                                                       12341, 12401, 25799,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47254, 0, 3,
                                                                       45994, 24859, 46084,
                                                                       12401, 12461, 25899,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47404, 0, 3,
                                                                       46084, 24919, 46174,
                                                                       12461, 12521, 25999,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 47554, 0, 3,
                                                                       46174, 24979, 46264,
                                                                       12521, 12581, 26099,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 47704, 0, 3,
                                                                       46354, 25299, 46504,
                                                                       12701, 12791, 26499,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 47929, 0, 3,
                                                                       46504, 25399, 46654,
                                                                       12791, 12881, 26649,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 48154, 0, 3,
                                                                       46654, 25499, 46804,
                                                                       12881, 12971, 26799,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 48379, 0, 3,
                                                                       46804, 25599, 46954,
                                                                       12971, 13061, 26949,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 48604, 0, 3,
                                                                       46954, 25699, 47104,
                                                                       13061, 13151, 27099,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 48829, 0, 3,
                                                                       47104, 25799, 47254,
                                                                       13151, 13241, 27249,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 49054, 0, 3,
                                                                       47254, 25899, 47404,
                                                                       13241, 13331, 27399,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 49279, 0, 3,
                                                                       47404, 25999, 47554,
                                                                       13331, 13421, 27549,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 49504, 0, 3,
                                                                       47704, 26499, 47929,
                                                                       13601, 13727, 28119,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 49819, 0, 3,
                                                                       47929, 26649, 48154,
                                                                       13727, 13853, 28329,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 50134, 0, 3,
                                                                       48154, 26799, 48379,
                                                                       13853, 13979, 28539,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 50449, 0, 3,
                                                                       48379, 26949, 48604,
                                                                       13979, 14105, 28749,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 50764, 0, 3,
                                                                       48604, 27099, 48829,
                                                                       14105, 14231, 28959,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 51079, 0, 3,
                                                                       48829, 27249, 49054,
                                                                       14231, 14357, 29169,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 51394, 0, 3,
                                                                       49054, 27399, 49279,
                                                                       14357, 14483, 29379,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 51709, 0, 3,
                                                                       49504, 28119, 49819,
                                                                       14735, 14903, 30149,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 52129, 0, 3,
                                                                       49819, 28329, 50134,
                                                                       14903, 15071, 30429,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 52549, 0, 3,
                                                                       50134, 28539, 50449,
                                                                       15071, 15239, 30709,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 52969, 0, 3,
                                                                       50449, 28749, 50764,
                                                                       15239, 15407, 30989,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 53389, 0, 3,
                                                                       50764, 28959, 51079,
                                                                       15407, 15575, 31269,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 53809, 0, 3,
                                                                       51079, 29169, 51394,
                                                                       15575, 15743, 31549,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 54229, 0, 3,
                                                                       51709, 30149, 52129,
                                                                       16079, 16295, 32549,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 54769, 0, 3,
                                                                       52129, 30429, 52549,
                                                                       16295, 16511, 32909,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 55309, 0, 3,
                                                                       52549, 30709, 52969,
                                                                       16511, 16727, 33269,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 55849, 0, 3,
                                                                       52969, 30989, 53389,
                                                                       16727, 16943, 33629,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 56389, 0, 3,
                                                                       53389, 31269, 53809,
                                                                       16943, 17159, 33989,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 56929, 0, 3,
                                                                       54229, 32549, 54769,
                                                                       17591, 17861, 35249,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 57604, 0, 3,
                                                                       54769, 32909, 55309,
                                                                       17861, 18131, 35699,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 58279, 0, 3,
                                                                       55309, 33269, 55849,
                                                                       18131, 18401, 36149,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 58954, 0, 3,
                                                                       55849, 33629, 56389,
                                                                       18401, 18671, 36599,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 59629, 0, 3,
                                                                       56929, 35249, 57604,
                                                                       19211, 19541, 38149,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 60454, 0, 3,
                                                                       57604, 35699, 58279,
                                                                       19541, 19871, 38699,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 61279, 0, 3,
                                                                       58279, 36149, 58954,
                                                                       19871, 20201, 39249,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 62104, 0, 3,
                                                                       59629, 38149, 60454,
                                                                       20861, 21257, 41119,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 63094, 0, 3,
                                                                       60454, 38699, 61279,
                                                                       21257, 21653, 41779,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 64084, 0, 3,
                                                                       62104, 41119, 63094,
                                                                       22445, 22913, 43999,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65254, 3, 23849,
                                                                       23859, 44779, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65275, 3, 23859,
                                                                       23869, 44794, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65296, 3, 23869,
                                                                       23879, 44809, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65317, 3, 23879,
                                                                       23889, 44824, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65338, 3, 23889,
                                                                       23899, 44839, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65359, 3, 23899,
                                                                       23909, 44854, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65380, 3, 23909,
                                                                       23919, 44869, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65401, 3, 23919,
                                                                       23929, 44884, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65422, 3, 23929,
                                                                       23939, 44899, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65443, 3, 23939,
                                                                       23949, 44914, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65464, 3, 23949,
                                                                       23959, 44929, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65485, 3, 23959,
                                                                       23969, 44944, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65506, 0, 3,
                                                                       65254, 44779, 65275,
                                                                       23989, 24019, 44959,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65569, 0, 3,
                                                                       65275, 44794, 65296,
                                                                       24019, 24049, 45004,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65632, 0, 3,
                                                                       65296, 44809, 65317,
                                                                       24049, 24079, 45049,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65695, 0, 3,
                                                                       65317, 44824, 65338,
                                                                       24079, 24109, 45094,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65758, 0, 3,
                                                                       65338, 44839, 65359,
                                                                       24109, 24139, 45139,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65821, 0, 3,
                                                                       65359, 44854, 65380,
                                                                       24139, 24169, 45184,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65884, 0, 3,
                                                                       65380, 44869, 65401,
                                                                       24169, 24199, 45229,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65947, 0, 3,
                                                                       65401, 44884, 65422,
                                                                       24199, 24229, 45274,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 66010, 0, 3,
                                                                       65422, 44899, 65443,
                                                                       24229, 24259, 45319,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 66073, 0, 3,
                                                                       65443, 44914, 65464,
                                                                       24259, 24289, 45364,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 66136, 0, 3,
                                                                       65464, 44929, 65485,
                                                                       24289, 24319, 45409,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 66199, 0, 3,
                                                                       65506, 44959, 65569,
                                                                       24379, 24439, 45454,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 66325, 0, 3,
                                                                       65569, 45004, 65632,
                                                                       24439, 24499, 45544,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 66451, 0, 3,
                                                                       65632, 45049, 65695,
                                                                       24499, 24559, 45634,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 66577, 0, 3,
                                                                       65695, 45094, 65758,
                                                                       24559, 24619, 45724,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 66703, 0, 3,
                                                                       65758, 45139, 65821,
                                                                       24619, 24679, 45814,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 66829, 0, 3,
                                                                       65821, 45184, 65884,
                                                                       24679, 24739, 45904,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 66955, 0, 3,
                                                                       65884, 45229, 65947,
                                                                       24739, 24799, 45994,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 67081, 0, 3,
                                                                       65947, 45274, 66010,
                                                                       24799, 24859, 46084,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 67207, 0, 3,
                                                                       66010, 45319, 66073,
                                                                       24859, 24919, 46174,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 67333, 0, 3,
                                                                       66073, 45364, 66136,
                                                                       24919, 24979, 46264,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 67459, 0, 3,
                                                                       66199, 45454, 66325,
                                                                       25099, 25199, 46354,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 67669, 0, 3,
                                                                       66325, 45544, 66451,
                                                                       25199, 25299, 46504,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 67879, 0, 3,
                                                                       66451, 45634, 66577,
                                                                       25299, 25399, 46654,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 68089, 0, 3,
                                                                       66577, 45724, 66703,
                                                                       25399, 25499, 46804,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 68299, 0, 3,
                                                                       66703, 45814, 66829,
                                                                       25499, 25599, 46954,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 68509, 0, 3,
                                                                       66829, 45904, 66955,
                                                                       25599, 25699, 47104,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 68719, 0, 3,
                                                                       66955, 45994, 67081,
                                                                       25699, 25799, 47254,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 68929, 0, 3,
                                                                       67081, 46084, 67207,
                                                                       25799, 25899, 47404,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 69139, 0, 3,
                                                                       67207, 46174, 67333,
                                                                       25899, 25999, 47554,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 69349, 0, 3,
                                                                       67459, 46354, 67669,
                                                                       26199, 26349, 47704,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 69664, 0, 3,
                                                                       67669, 46504, 67879,
                                                                       26349, 26499, 47929,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 69979, 0, 3,
                                                                       67879, 46654, 68089,
                                                                       26499, 26649, 48154,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 70294, 0, 3,
                                                                       68089, 46804, 68299,
                                                                       26649, 26799, 48379,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 70609, 0, 3,
                                                                       68299, 46954, 68509,
                                                                       26799, 26949, 48604,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 70924, 0, 3,
                                                                       68509, 47104, 68719,
                                                                       26949, 27099, 48829,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 71239, 0, 3,
                                                                       68719, 47254, 68929,
                                                                       27099, 27249, 49054,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 71554, 0, 3,
                                                                       68929, 47404, 69139,
                                                                       27249, 27399, 49279,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 71869, 0, 3,
                                                                       69349, 47704, 69664,
                                                                       27699, 27909, 49504,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 72310, 0, 3,
                                                                       69664, 47929, 69979,
                                                                       27909, 28119, 49819,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 72751, 0, 3,
                                                                       69979, 48154, 70294,
                                                                       28119, 28329, 50134,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 73192, 0, 3,
                                                                       70294, 48379, 70609,
                                                                       28329, 28539, 50449,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 73633, 0, 3,
                                                                       70609, 48604, 70924,
                                                                       28539, 28749, 50764,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 74074, 0, 3,
                                                                       70924, 48829, 71239,
                                                                       28749, 28959, 51079,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 74515, 0, 3,
                                                                       71239, 49054, 71554,
                                                                       28959, 29169, 51394,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 74956, 0, 3,
                                                                       71869, 49504, 72310,
                                                                       29589, 29869, 51709,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 75544, 0, 3,
                                                                       72310, 49819, 72751,
                                                                       29869, 30149, 52129,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 76132, 0, 3,
                                                                       72751, 50134, 73192,
                                                                       30149, 30429, 52549,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 76720, 0, 3,
                                                                       73192, 50449, 73633,
                                                                       30429, 30709, 52969,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 77308, 0, 3,
                                                                       73633, 50764, 74074,
                                                                       30709, 30989, 53389,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 77896, 0, 3,
                                                                       74074, 51079, 74515,
                                                                       30989, 31269, 53809,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 78484, 0, 3,
                                                                       74956, 51709, 75544,
                                                                       31829, 32189, 54229,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 79240, 0, 3,
                                                                       75544, 52129, 76132,
                                                                       32189, 32549, 54769,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 79996, 0, 3,
                                                                       76132, 52549, 76720,
                                                                       32549, 32909, 55309,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 80752, 0, 3,
                                                                       76720, 52969, 77308,
                                                                       32909, 33269, 55849,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 81508, 0, 3,
                                                                       77308, 53389, 77896,
                                                                       33269, 33629, 56389,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 82264, 0, 3,
                                                                       78484, 54229, 79240,
                                                                       34349, 34799, 56929,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 83209, 0, 3,
                                                                       79240, 54769, 79996,
                                                                       34799, 35249, 57604,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 84154, 0, 3,
                                                                       79996, 55309, 80752,
                                                                       35249, 35699, 58279,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 85099, 0, 3,
                                                                       80752, 55849, 81508,
                                                                       35699, 36149, 58954,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 86044, 0, 3,
                                                                       82264, 56929, 83209,
                                                                       37049, 37599, 59629,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 87199, 0, 3,
                                                                       83209, 57604, 84154,
                                                                       37599, 38149, 60454,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 88354, 0, 3,
                                                                       84154, 58279, 85099,
                                                                       38149, 38699, 61279,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 89509, 0, 3,
                                                                       86044, 59629, 87199,
                                                                       39799, 40459, 62104,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 90895, 0, 3,
                                                                       87199, 60454, 88354,
                                                                       40459, 41119, 63094,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 92281, 0, 3,
                                                                       89509, 62104, 90895,
                                                                       42439, 43219, 64084,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 93919, 74956, 588, ncols);

                    simdfunc::contract_primitives(buffer, 94815, 78484, 756, ncols);

                    simdfunc::contract_primitives(buffer, 95967, 82264, 945, ncols);

                    simdfunc::contract_primitives(buffer, 97407, 86044, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 99167, 89509, 1386, ncols);

                    simdfunc::contract_primitives(buffer, 101279, 92281, 1638, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 94507, 93919, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 95571, 94815, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 96912, 95967, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 98562, 97407, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 100553, 99167, 66, 1, nmax);

        simdtrf::transform_h_inner(buffer, 102917, 101279, 78, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 103775, 94507, 95571, 11,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 104699, 95571, 96912, 11,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 105887, 96912, 98562, 11,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 107372, 98562, 100553, 11,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 109187, 100553, 102917, 11,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 111365, 103775, 104699, 11,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 113213, 104699, 105887, 11,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 115589, 105887, 107372, 11,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 118559, 107372, 109187, 11,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 122189, 111365, 113213, 11,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 125269, 113213, 115589, 11,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 129229, 115589, 118559, 11,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 134179, 122189, 125269, 11,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 138799, 125269, 129229, 11,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 144739, 134179, 138799, 11,
                                             nmax);

        simdtrf::transform_h_inner(buffer, 151207, 144739, 28, 11, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 151207, 121, nmax);
    }

    for (size_t m = 0; m < 1573; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
