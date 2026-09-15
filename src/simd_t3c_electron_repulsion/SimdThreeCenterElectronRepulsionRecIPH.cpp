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


#include "SimdThreeCenterElectronRepulsionRecIPH.hpp"

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
#include "SimdTransferIP.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_iph_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_iph_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 28144, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 429 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 28144, 24248, 1652, dimensions);

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

                simdfunc::compute_pair_exponent(buffer, coordinates, 6, nmax, mu);

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 7, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12}, ncols, fj, 6,
                                                        fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 20, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 23, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 26, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 29, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 32, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 35, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 38, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 41, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 44, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 47, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 50, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 53, 0, 3, 8, 9,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 59, 0, 3, 9, 10,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 65, 0, 3, 10, 11,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 71, 0, 3, 11, 12,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 77, 0, 3, 12, 13,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 83, 0, 3, 13, 14,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 89, 0, 3, 14, 15,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 95, 0, 3, 15, 16,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 101, 0, 3, 16, 17,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 107, 0, 3, 17, 18,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 113, 0, 3, 20, 23,
                                                                       53, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 123, 0, 3, 23, 26,
                                                                       59, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 133, 0, 3, 26, 29,
                                                                       65, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 143, 0, 3, 29, 32,
                                                                       71, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 153, 0, 3, 32, 35,
                                                                       77, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 163, 0, 3, 35, 38,
                                                                       83, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 173, 0, 3, 38, 41,
                                                                       89, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 183, 0, 3, 41, 44,
                                                                       95, 101, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 193, 0, 3, 44, 47,
                                                                       101, 107, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 203, 0, 3, 53, 59,
                                                                       113, 123, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 218, 0, 3, 59, 65,
                                                                       123, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 233, 0, 3, 65, 71,
                                                                       133, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 248, 0, 3, 71, 77,
                                                                       143, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 263, 0, 3, 77, 83,
                                                                       153, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 278, 0, 3, 83, 89,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 293, 0, 3, 89, 95,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 308, 0, 3, 95,
                                                                       101, 183, 193, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 323, 0, 3, 113,
                                                                       123, 203, 218, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 344, 0, 3, 123,
                                                                       133, 218, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 365, 0, 3, 133,
                                                                       143, 233, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 386, 0, 3, 143,
                                                                       153, 248, 263, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 407, 0, 3, 153,
                                                                       163, 263, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 163,
                                                                       173, 278, 293, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 449, 0, 3, 173,
                                                                       183, 293, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 470, 0, 3, 203,
                                                                       218, 323, 344, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 498, 0, 3, 218,
                                                                       233, 344, 365, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 526, 0, 3, 233,
                                                                       248, 365, 386, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 554, 0, 3, 248,
                                                                       263, 386, 407, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 582, 0, 3, 263,
                                                                       278, 407, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 610, 0, 3, 278,
                                                                       293, 428, 449, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 638, 0, 3, 323,
                                                                       344, 470, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 674, 0, 3, 344,
                                                                       365, 498, 526, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 710, 0, 3, 365,
                                                                       386, 526, 554, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 746, 0, 3, 386,
                                                                       407, 554, 582, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 782, 0, 3, 407,
                                                                       428, 582, 610, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 818, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 821, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 824, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 827, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 830, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 833, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 836, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 839, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 842, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 845, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 848, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 851, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 854, 3, 8, 20,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 863, 3, 9, 23,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 872, 3, 10, 26,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 881, 3, 11, 29,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 890, 3, 12, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 899, 3, 13, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 908, 3, 14, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 917, 3, 15, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 926, 3, 16, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 935, 3, 17, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 944, 3, 18, 50,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 953, 3, 20, 53,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 971, 3, 23, 59,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 989, 3, 26, 65,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1007, 3, 29, 71,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1025, 3, 32, 77,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1043, 3, 35, 83,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1061, 3, 38, 89,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1079, 3, 41, 95,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1097, 3, 44, 101,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1115, 3, 47, 107,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1133, 3, 53, 113,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1163, 3, 59, 123,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1193, 3, 65, 133,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1223, 3, 71, 143,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1253, 3, 77, 153,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1283, 3, 83, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1313, 3, 89, 173,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1343, 3, 95, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1373, 3, 101, 193,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1403, 3, 113, 203,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1448, 3, 123, 218,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1493, 3, 133, 233,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1538, 3, 143, 248,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1583, 3, 153, 263,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1628, 3, 163, 278,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1673, 3, 173, 293,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1718, 3, 183, 308,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1763, 3, 203, 323,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1826, 3, 218, 344,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1889, 3, 233, 365,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1952, 3, 248, 386,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2015, 3, 263, 407,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2078, 3, 278, 428,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2141, 3, 293, 449,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2204, 3, 323, 470,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2288, 3, 344, 498,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2372, 3, 365, 526,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2456, 3, 386, 554,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2540, 3, 407, 582,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2624, 3, 428, 610,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 2708, 3, 470, 638,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 2816, 3, 498, 674,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 2924, 3, 526, 710,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3032, 3, 554, 746,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3140, 3, 582, 782,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3248, 3, 8, 9,
                                                                       824, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3254, 3, 9, 10,
                                                                       827, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3260, 3, 10, 11,
                                                                       830, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3266, 3, 11, 12,
                                                                       833, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3272, 3, 12, 13,
                                                                       836, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3278, 3, 13, 14,
                                                                       839, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3284, 3, 14, 15,
                                                                       842, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3290, 3, 15, 16,
                                                                       845, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3296, 3, 16, 17,
                                                                       848, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3302, 3, 17, 18,
                                                                       851, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3308, 0, 3, 3248,
                                                                       824, 3254, 872, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3326, 0, 3, 3254,
                                                                       827, 3260, 881, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3344, 0, 3, 3260,
                                                                       830, 3266, 890, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3362, 0, 3, 3266,
                                                                       833, 3272, 899, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3380, 0, 3, 3272,
                                                                       836, 3278, 908, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3398, 0, 3, 3278,
                                                                       839, 3284, 917, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3416, 0, 3, 3284,
                                                                       842, 3290, 926, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3434, 0, 3, 3290,
                                                                       845, 3296, 935, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3452, 0, 3, 3296,
                                                                       848, 3302, 944, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3470, 0, 3, 3308,
                                                                       872, 3326, 53, 59, 989,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3506, 0, 3, 3326,
                                                                       881, 3344, 59, 65, 1007,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3542, 0, 3, 3344,
                                                                       890, 3362, 65, 71, 1025,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3578, 0, 3, 3362,
                                                                       899, 3380, 71, 77, 1043,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3614, 0, 3, 3380,
                                                                       908, 3398, 77, 83, 1061,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3650, 0, 3, 3398,
                                                                       917, 3416, 83, 89, 1079,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3686, 0, 3, 3416,
                                                                       926, 3434, 89, 95, 1097,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3722, 0, 3, 3434,
                                                                       935, 3452, 95, 101, 1115,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3758, 0, 3, 3470,
                                                                       989, 3506, 113, 123, 1193,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3818, 0, 3, 3506,
                                                                       1007, 3542, 123, 133,
                                                                       1223, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3878, 0, 3, 3542,
                                                                       1025, 3578, 133, 143,
                                                                       1253, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3938, 0, 3, 3578,
                                                                       1043, 3614, 143, 153,
                                                                       1283, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3998, 0, 3, 3614,
                                                                       1061, 3650, 153, 163,
                                                                       1313, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4058, 0, 3, 3650,
                                                                       1079, 3686, 163, 173,
                                                                       1343, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4118, 0, 3, 3686,
                                                                       1097, 3722, 173, 183,
                                                                       1373, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4178, 0, 3, 3758,
                                                                       1193, 3818, 203, 218,
                                                                       1493, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4268, 0, 3, 3818,
                                                                       1223, 3878, 218, 233,
                                                                       1538, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4358, 0, 3, 3878,
                                                                       1253, 3938, 233, 248,
                                                                       1583, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4448, 0, 3, 3938,
                                                                       1283, 3998, 248, 263,
                                                                       1628, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4538, 0, 3, 3998,
                                                                       1313, 4058, 263, 278,
                                                                       1673, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4628, 0, 3, 4058,
                                                                       1343, 4118, 278, 293,
                                                                       1718, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4718, 0, 3, 4178,
                                                                       1493, 4268, 323, 344,
                                                                       1889, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4844, 0, 3, 4268,
                                                                       1538, 4358, 344, 365,
                                                                       1952, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 4970, 0, 3, 4358,
                                                                       1583, 4448, 365, 386,
                                                                       2015, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5096, 0, 3, 4448,
                                                                       1628, 4538, 386, 407,
                                                                       2078, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5222, 0, 3, 4538,
                                                                       1673, 4628, 407, 428,
                                                                       2141, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 5348, 0, 3, 4718,
                                                                       1889, 4844, 470, 498,
                                                                       2372, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 5516, 0, 3, 4844,
                                                                       1952, 4970, 498, 526,
                                                                       2456, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 5684, 0, 3, 4970,
                                                                       2015, 5096, 526, 554,
                                                                       2540, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 5852, 0, 3, 5096,
                                                                       2078, 5222, 554, 582,
                                                                       2624, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 6020, 0, 3, 5348,
                                                                       2372, 5516, 638, 674,
                                                                       2924, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 6236, 0, 3, 5516,
                                                                       2456, 5684, 674, 710,
                                                                       3032, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 6452, 0, 3, 5684,
                                                                       2540, 5852, 710, 746,
                                                                       3140, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6668, 3, 818, 821,
                                                                       3248, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6678, 3, 821, 824,
                                                                       3254, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6688, 3, 824, 827,
                                                                       3260, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6698, 3, 827, 830,
                                                                       3266, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6708, 3, 830, 833,
                                                                       3272, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6718, 3, 833, 836,
                                                                       3278, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6728, 3, 836, 839,
                                                                       3284, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6738, 3, 839, 842,
                                                                       3290, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6748, 3, 842, 845,
                                                                       3296, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6758, 3, 845, 848,
                                                                       3302, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6768, 0, 3, 6668,
                                                                       3248, 6678, 854, 863,
                                                                       3308, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6798, 0, 3, 6678,
                                                                       3254, 6688, 863, 872,
                                                                       3326, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6828, 0, 3, 6688,
                                                                       3260, 6698, 872, 881,
                                                                       3344, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6858, 0, 3, 6698,
                                                                       3266, 6708, 881, 890,
                                                                       3362, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6888, 0, 3, 6708,
                                                                       3272, 6718, 890, 899,
                                                                       3380, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6918, 0, 3, 6718,
                                                                       3278, 6728, 899, 908,
                                                                       3398, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6948, 0, 3, 6728,
                                                                       3284, 6738, 908, 917,
                                                                       3416, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6978, 0, 3, 6738,
                                                                       3290, 6748, 917, 926,
                                                                       3434, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7008, 0, 3, 6748,
                                                                       3296, 6758, 926, 935,
                                                                       3452, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7038, 0, 3, 6768,
                                                                       3308, 6798, 953, 971,
                                                                       3470, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7098, 0, 3, 6798,
                                                                       3326, 6828, 971, 989,
                                                                       3506, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7158, 0, 3, 6828,
                                                                       3344, 6858, 989, 1007,
                                                                       3542, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7218, 0, 3, 6858,
                                                                       3362, 6888, 1007, 1025,
                                                                       3578, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7278, 0, 3, 6888,
                                                                       3380, 6918, 1025, 1043,
                                                                       3614, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7338, 0, 3, 6918,
                                                                       3398, 6948, 1043, 1061,
                                                                       3650, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7398, 0, 3, 6948,
                                                                       3416, 6978, 1061, 1079,
                                                                       3686, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7458, 0, 3, 6978,
                                                                       3434, 7008, 1079, 1097,
                                                                       3722, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7518, 0, 3, 7038,
                                                                       3470, 7098, 1133, 1163,
                                                                       3758, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7618, 0, 3, 7098,
                                                                       3506, 7158, 1163, 1193,
                                                                       3818, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7718, 0, 3, 7158,
                                                                       3542, 7218, 1193, 1223,
                                                                       3878, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7818, 0, 3, 7218,
                                                                       3578, 7278, 1223, 1253,
                                                                       3938, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7918, 0, 3, 7278,
                                                                       3614, 7338, 1253, 1283,
                                                                       3998, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 8018, 0, 3, 7338,
                                                                       3650, 7398, 1283, 1313,
                                                                       4058, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 8118, 0, 3, 7398,
                                                                       3686, 7458, 1313, 1343,
                                                                       4118, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8218, 0, 3, 7518,
                                                                       3758, 7618, 1403, 1448,
                                                                       4178, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8368, 0, 3, 7618,
                                                                       3818, 7718, 1448, 1493,
                                                                       4268, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8518, 0, 3, 7718,
                                                                       3878, 7818, 1493, 1538,
                                                                       4358, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8668, 0, 3, 7818,
                                                                       3938, 7918, 1538, 1583,
                                                                       4448, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8818, 0, 3, 7918,
                                                                       3998, 8018, 1583, 1628,
                                                                       4538, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8968, 0, 3, 8018,
                                                                       4058, 8118, 1628, 1673,
                                                                       4628, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 9118, 0, 3, 8218,
                                                                       4178, 8368, 1763, 1826,
                                                                       4718, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 9328, 0, 3, 8368,
                                                                       4268, 8518, 1826, 1889,
                                                                       4844, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 9538, 0, 3, 8518,
                                                                       4358, 8668, 1889, 1952,
                                                                       4970, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 9748, 0, 3, 8668,
                                                                       4448, 8818, 1952, 2015,
                                                                       5096, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 9958, 0, 3, 8818,
                                                                       4538, 8968, 2015, 2078,
                                                                       5222, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 10168, 0, 3, 9118,
                                                                       4718, 9328, 2204, 2288,
                                                                       5348, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 10448, 0, 3, 9328,
                                                                       4844, 9538, 2288, 2372,
                                                                       5516, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 10728, 0, 3, 9538,
                                                                       4970, 9748, 2372, 2456,
                                                                       5684, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 11008, 0, 3, 9748,
                                                                       5096, 9958, 2456, 2540,
                                                                       5852, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 11288, 0, 3,
                                                                       10168, 5348, 10448, 2708,
                                                                       2816, 6020, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 11648, 0, 3,
                                                                       10448, 5516, 10728, 2816,
                                                                       2924, 6236, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 12008, 0, 3,
                                                                       10728, 5684, 11008, 2924,
                                                                       3032, 6452, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12368, 3, 3248,
                                                                       3254, 6688, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12383, 3, 3254,
                                                                       3260, 6698, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12398, 3, 3260,
                                                                       3266, 6708, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12413, 3, 3266,
                                                                       3272, 6718, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12428, 3, 3272,
                                                                       3278, 6728, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12443, 3, 3278,
                                                                       3284, 6738, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12458, 3, 3284,
                                                                       3290, 6748, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12473, 3, 3290,
                                                                       3296, 6758, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 12488, 0, 3,
                                                                       12368, 6688, 12383, 3308,
                                                                       3326, 6828, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 12533, 0, 3,
                                                                       12383, 6698, 12398, 3326,
                                                                       3344, 6858, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 12578, 0, 3,
                                                                       12398, 6708, 12413, 3344,
                                                                       3362, 6888, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 12623, 0, 3,
                                                                       12413, 6718, 12428, 3362,
                                                                       3380, 6918, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 12668, 0, 3,
                                                                       12428, 6728, 12443, 3380,
                                                                       3398, 6948, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 12713, 0, 3,
                                                                       12443, 6738, 12458, 3398,
                                                                       3416, 6978, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 12758, 0, 3,
                                                                       12458, 6748, 12473, 3416,
                                                                       3434, 7008, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 12803, 0, 3,
                                                                       12488, 6828, 12533, 3470,
                                                                       3506, 7158, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 12893, 0, 3,
                                                                       12533, 6858, 12578, 3506,
                                                                       3542, 7218, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 12983, 0, 3,
                                                                       12578, 6888, 12623, 3542,
                                                                       3578, 7278, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 13073, 0, 3,
                                                                       12623, 6918, 12668, 3578,
                                                                       3614, 7338, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 13163, 0, 3,
                                                                       12668, 6948, 12713, 3614,
                                                                       3650, 7398, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 13253, 0, 3,
                                                                       12713, 6978, 12758, 3650,
                                                                       3686, 7458, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13343, 0, 3,
                                                                       12803, 7158, 12893, 3758,
                                                                       3818, 7718, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13493, 0, 3,
                                                                       12893, 7218, 12983, 3818,
                                                                       3878, 7818, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13643, 0, 3,
                                                                       12983, 7278, 13073, 3878,
                                                                       3938, 7918, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13793, 0, 3,
                                                                       13073, 7338, 13163, 3938,
                                                                       3998, 8018, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13943, 0, 3,
                                                                       13163, 7398, 13253, 3998,
                                                                       4058, 8118, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 14093, 0, 3,
                                                                       13343, 7718, 13493, 4178,
                                                                       4268, 8518, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 14318, 0, 3,
                                                                       13493, 7818, 13643, 4268,
                                                                       4358, 8668, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 14543, 0, 3,
                                                                       13643, 7918, 13793, 4358,
                                                                       4448, 8818, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 14768, 0, 3,
                                                                       13793, 8018, 13943, 4448,
                                                                       4538, 8968, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 14993, 0, 3,
                                                                       14093, 8518, 14318, 4718,
                                                                       4844, 9538, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 15308, 0, 3,
                                                                       14318, 8668, 14543, 4844,
                                                                       4970, 9748, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 15623, 0, 3,
                                                                       14543, 8818, 14768, 4970,
                                                                       5096, 9958, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 15938, 0, 3,
                                                                       14993, 9538, 15308, 5348,
                                                                       5516, 10728, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 16358, 0, 3,
                                                                       15308, 9748, 15623, 5516,
                                                                       5684, 11008, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 16778, 0, 3,
                                                                       15938, 10728, 16358, 6020,
                                                                       6236, 12008, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17318, 3, 6668,
                                                                       6678, 12368, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17339, 3, 6678,
                                                                       6688, 12383, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17360, 3, 6688,
                                                                       6698, 12398, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17381, 3, 6698,
                                                                       6708, 12413, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17402, 3, 6708,
                                                                       6718, 12428, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17423, 3, 6718,
                                                                       6728, 12443, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17444, 3, 6728,
                                                                       6738, 12458, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17465, 3, 6738,
                                                                       6748, 12473, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17486, 0, 3,
                                                                       17318, 12368, 17339, 6768,
                                                                       6798, 12488, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17549, 0, 3,
                                                                       17339, 12383, 17360, 6798,
                                                                       6828, 12533, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17612, 0, 3,
                                                                       17360, 12398, 17381, 6828,
                                                                       6858, 12578, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17675, 0, 3,
                                                                       17381, 12413, 17402, 6858,
                                                                       6888, 12623, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17738, 0, 3,
                                                                       17402, 12428, 17423, 6888,
                                                                       6918, 12668, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17801, 0, 3,
                                                                       17423, 12443, 17444, 6918,
                                                                       6948, 12713, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17864, 0, 3,
                                                                       17444, 12458, 17465, 6948,
                                                                       6978, 12758, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 17927, 0, 3,
                                                                       17486, 12488, 17549, 7038,
                                                                       7098, 12803, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18053, 0, 3,
                                                                       17549, 12533, 17612, 7098,
                                                                       7158, 12893, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18179, 0, 3,
                                                                       17612, 12578, 17675, 7158,
                                                                       7218, 12983, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18305, 0, 3,
                                                                       17675, 12623, 17738, 7218,
                                                                       7278, 13073, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18431, 0, 3,
                                                                       17738, 12668, 17801, 7278,
                                                                       7338, 13163, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18557, 0, 3,
                                                                       17801, 12713, 17864, 7338,
                                                                       7398, 13253, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 18683, 0, 3,
                                                                       17927, 12803, 18053, 7518,
                                                                       7618, 13343, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 18893, 0, 3,
                                                                       18053, 12893, 18179, 7618,
                                                                       7718, 13493, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 19103, 0, 3,
                                                                       18179, 12983, 18305, 7718,
                                                                       7818, 13643, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 19313, 0, 3,
                                                                       18305, 13073, 18431, 7818,
                                                                       7918, 13793, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 19523, 0, 3,
                                                                       18431, 13163, 18557, 7918,
                                                                       8018, 13943, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 19733, 0, 3,
                                                                       18683, 13343, 18893, 8218,
                                                                       8368, 14093, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 20048, 0, 3,
                                                                       18893, 13493, 19103, 8368,
                                                                       8518, 14318, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 20363, 0, 3,
                                                                       19103, 13643, 19313, 8518,
                                                                       8668, 14543, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 20678, 0, 3,
                                                                       19313, 13793, 19523, 8668,
                                                                       8818, 14768, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 20993, 0, 3,
                                                                       19733, 14093, 20048, 9118,
                                                                       9328, 14993, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 21434, 0, 3,
                                                                       20048, 14318, 20363, 9328,
                                                                       9538, 15308, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 21875, 0, 3,
                                                                       20363, 14543, 20678, 9538,
                                                                       9748, 15623, ncols, gamma,
                                                                       p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 22316, 0, 3,
                                                                       20993, 14993, 21434,
                                                                       10168, 10448, 15938,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 22904, 0, 3,
                                                                       21434, 15308, 21875,
                                                                       10448, 10728, 16358,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 23492, 0, 3,
                                                                       22316, 15938, 22904,
                                                                       11288, 11648, 16778,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 24248, 22316, 588, ncols);

                    simdfunc::contract_primitives(buffer, 25144, 23492, 756, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 24836, 24248, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 25900, 25144, 36, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 26296, 24836, 25900, 11,
                                             nmax);

        simdtrf::transform_p_inner(buffer, 27220, 26296, 28, 11, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 27220, 33, nmax);
    }

    for (size_t m = 0; m < 429; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
