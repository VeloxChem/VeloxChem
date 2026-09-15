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


#include "SimdThreeCenterElectronRepulsionRecFDK.hpp"

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
#include "SimdTransferFD.hpp"
#include "SimdTransferFP.hpp"
#include "SimdTransferGP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_fdk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_fdk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 29209, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 525 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 29209, 24088, 2031, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 470, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 473, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 476, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 479, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 482, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 485, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 488, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 491, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 494, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 497, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 500, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 503, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 506, 3, 8, 20,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 515, 3, 9, 23,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 524, 3, 10, 26,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 533, 3, 11, 29,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 542, 3, 12, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 551, 3, 13, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 560, 3, 14, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 569, 3, 15, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 578, 3, 16, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 587, 3, 17, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 596, 3, 18, 50,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 605, 3, 20, 53,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 623, 3, 23, 59,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 641, 3, 26, 65,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 659, 3, 29, 71,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 677, 3, 32, 77,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 695, 3, 35, 83,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 713, 3, 38, 89,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 731, 3, 41, 95,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 749, 3, 44, 101,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 767, 3, 47, 107,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 785, 3, 53, 113,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 815, 3, 59, 123,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 845, 3, 65, 133,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 875, 3, 71, 143,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 905, 3, 77, 153,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 935, 3, 83, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 965, 3, 89, 173,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 995, 3, 95, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1025, 3, 101, 193,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1055, 3, 113, 203,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1100, 3, 123, 218,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1145, 3, 133, 233,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1190, 3, 143, 248,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1235, 3, 153, 263,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1280, 3, 163, 278,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1325, 3, 173, 293,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1370, 3, 183, 308,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1415, 3, 203, 323,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1478, 3, 218, 344,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1541, 3, 233, 365,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1604, 3, 248, 386,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1667, 3, 263, 407,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1730, 3, 278, 428,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1793, 3, 293, 449,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1856, 3, 8, 9,
                                                                       476, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1862, 3, 9, 10,
                                                                       479, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1868, 3, 10, 11,
                                                                       482, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1874, 3, 11, 12,
                                                                       485, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1880, 3, 12, 13,
                                                                       488, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1886, 3, 13, 14,
                                                                       491, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1892, 3, 14, 15,
                                                                       494, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1898, 3, 15, 16,
                                                                       497, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1904, 3, 16, 17,
                                                                       500, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1910, 3, 17, 18,
                                                                       503, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1916, 0, 3, 1856,
                                                                       476, 1862, 524, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1934, 0, 3, 1862,
                                                                       479, 1868, 533, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1952, 0, 3, 1868,
                                                                       482, 1874, 542, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1970, 0, 3, 1874,
                                                                       485, 1880, 551, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1988, 0, 3, 1880,
                                                                       488, 1886, 560, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2006, 0, 3, 1886,
                                                                       491, 1892, 569, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2024, 0, 3, 1892,
                                                                       494, 1898, 578, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2042, 0, 3, 1898,
                                                                       497, 1904, 587, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2060, 0, 3, 1904,
                                                                       500, 1910, 596, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2078, 0, 3, 1916,
                                                                       524, 1934, 53, 59, 641,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2114, 0, 3, 1934,
                                                                       533, 1952, 59, 65, 659,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2150, 0, 3, 1952,
                                                                       542, 1970, 65, 71, 677,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2186, 0, 3, 1970,
                                                                       551, 1988, 71, 77, 695,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2222, 0, 3, 1988,
                                                                       560, 2006, 77, 83, 713,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2258, 0, 3, 2006,
                                                                       569, 2024, 83, 89, 731,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2294, 0, 3, 2024,
                                                                       578, 2042, 89, 95, 749,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2330, 0, 3, 2042,
                                                                       587, 2060, 95, 101, 767,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2366, 0, 3, 2078,
                                                                       641, 2114, 113, 123, 845,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2426, 0, 3, 2114,
                                                                       659, 2150, 123, 133, 875,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2486, 0, 3, 2150,
                                                                       677, 2186, 133, 143, 905,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2546, 0, 3, 2186,
                                                                       695, 2222, 143, 153, 935,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2606, 0, 3, 2222,
                                                                       713, 2258, 153, 163, 965,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2666, 0, 3, 2258,
                                                                       731, 2294, 163, 173, 995,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2726, 0, 3, 2294,
                                                                       749, 2330, 173, 183, 1025,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2786, 0, 3, 2366,
                                                                       845, 2426, 203, 218, 1145,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2876, 0, 3, 2426,
                                                                       875, 2486, 218, 233, 1190,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2966, 0, 3, 2486,
                                                                       905, 2546, 233, 248, 1235,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3056, 0, 3, 2546,
                                                                       935, 2606, 248, 263, 1280,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3146, 0, 3, 2606,
                                                                       965, 2666, 263, 278, 1325,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3236, 0, 3, 2666,
                                                                       995, 2726, 278, 293, 1370,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 3326, 0, 3, 2786,
                                                                       1145, 2876, 323, 344,
                                                                       1541, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 3452, 0, 3, 2876,
                                                                       1190, 2966, 344, 365,
                                                                       1604, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 3578, 0, 3, 2966,
                                                                       1235, 3056, 365, 386,
                                                                       1667, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 3704, 0, 3, 3056,
                                                                       1280, 3146, 386, 407,
                                                                       1730, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 3830, 0, 3, 3146,
                                                                       1325, 3236, 407, 428,
                                                                       1793, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3956, 3, 470, 473,
                                                                       1856, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3966, 3, 473, 476,
                                                                       1862, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3976, 3, 476, 479,
                                                                       1868, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3986, 3, 479, 482,
                                                                       1874, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3996, 3, 482, 485,
                                                                       1880, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4006, 3, 485, 488,
                                                                       1886, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4016, 3, 488, 491,
                                                                       1892, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4026, 3, 491, 494,
                                                                       1898, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4036, 3, 494, 497,
                                                                       1904, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4046, 3, 497, 500,
                                                                       1910, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4056, 0, 3, 3956,
                                                                       1856, 3966, 506, 515,
                                                                       1916, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4086, 0, 3, 3966,
                                                                       1862, 3976, 515, 524,
                                                                       1934, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4116, 0, 3, 3976,
                                                                       1868, 3986, 524, 533,
                                                                       1952, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4146, 0, 3, 3986,
                                                                       1874, 3996, 533, 542,
                                                                       1970, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4176, 0, 3, 3996,
                                                                       1880, 4006, 542, 551,
                                                                       1988, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4206, 0, 3, 4006,
                                                                       1886, 4016, 551, 560,
                                                                       2006, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4236, 0, 3, 4016,
                                                                       1892, 4026, 560, 569,
                                                                       2024, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4266, 0, 3, 4026,
                                                                       1898, 4036, 569, 578,
                                                                       2042, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4296, 0, 3, 4036,
                                                                       1904, 4046, 578, 587,
                                                                       2060, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4326, 0, 3, 4056,
                                                                       1916, 4086, 605, 623,
                                                                       2078, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4386, 0, 3, 4086,
                                                                       1934, 4116, 623, 641,
                                                                       2114, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4446, 0, 3, 4116,
                                                                       1952, 4146, 641, 659,
                                                                       2150, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4506, 0, 3, 4146,
                                                                       1970, 4176, 659, 677,
                                                                       2186, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4566, 0, 3, 4176,
                                                                       1988, 4206, 677, 695,
                                                                       2222, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4626, 0, 3, 4206,
                                                                       2006, 4236, 695, 713,
                                                                       2258, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4686, 0, 3, 4236,
                                                                       2024, 4266, 713, 731,
                                                                       2294, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4746, 0, 3, 4266,
                                                                       2042, 4296, 731, 749,
                                                                       2330, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4806, 0, 3, 4326,
                                                                       2078, 4386, 785, 815,
                                                                       2366, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4906, 0, 3, 4386,
                                                                       2114, 4446, 815, 845,
                                                                       2426, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 5006, 0, 3, 4446,
                                                                       2150, 4506, 845, 875,
                                                                       2486, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 5106, 0, 3, 4506,
                                                                       2186, 4566, 875, 905,
                                                                       2546, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 5206, 0, 3, 4566,
                                                                       2222, 4626, 905, 935,
                                                                       2606, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 5306, 0, 3, 4626,
                                                                       2258, 4686, 935, 965,
                                                                       2666, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 5406, 0, 3, 4686,
                                                                       2294, 4746, 965, 995,
                                                                       2726, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 5506, 0, 3, 4806,
                                                                       2366, 4906, 1055, 1100,
                                                                       2786, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 5656, 0, 3, 4906,
                                                                       2426, 5006, 1100, 1145,
                                                                       2876, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 5806, 0, 3, 5006,
                                                                       2486, 5106, 1145, 1190,
                                                                       2966, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 5956, 0, 3, 5106,
                                                                       2546, 5206, 1190, 1235,
                                                                       3056, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 6106, 0, 3, 5206,
                                                                       2606, 5306, 1235, 1280,
                                                                       3146, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 6256, 0, 3, 5306,
                                                                       2666, 5406, 1280, 1325,
                                                                       3236, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 6406, 0, 3, 5506,
                                                                       2786, 5656, 1415, 1478,
                                                                       3326, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 6616, 0, 3, 5656,
                                                                       2876, 5806, 1478, 1541,
                                                                       3452, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 6826, 0, 3, 5806,
                                                                       2966, 5956, 1541, 1604,
                                                                       3578, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 7036, 0, 3, 5956,
                                                                       3056, 6106, 1604, 1667,
                                                                       3704, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 7246, 0, 3, 6106,
                                                                       3146, 6256, 1667, 1730,
                                                                       3830, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7456, 3, 1856,
                                                                       1862, 3976, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7471, 3, 1862,
                                                                       1868, 3986, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7486, 3, 1868,
                                                                       1874, 3996, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7501, 3, 1874,
                                                                       1880, 4006, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7516, 3, 1880,
                                                                       1886, 4016, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7531, 3, 1886,
                                                                       1892, 4026, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7546, 3, 1892,
                                                                       1898, 4036, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7561, 3, 1898,
                                                                       1904, 4046, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7576, 0, 3, 7456,
                                                                       3976, 7471, 1916, 1934,
                                                                       4116, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7621, 0, 3, 7471,
                                                                       3986, 7486, 1934, 1952,
                                                                       4146, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7666, 0, 3, 7486,
                                                                       3996, 7501, 1952, 1970,
                                                                       4176, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7711, 0, 3, 7501,
                                                                       4006, 7516, 1970, 1988,
                                                                       4206, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7756, 0, 3, 7516,
                                                                       4016, 7531, 1988, 2006,
                                                                       4236, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7801, 0, 3, 7531,
                                                                       4026, 7546, 2006, 2024,
                                                                       4266, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7846, 0, 3, 7546,
                                                                       4036, 7561, 2024, 2042,
                                                                       4296, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 7891, 0, 3, 7576,
                                                                       4116, 7621, 2078, 2114,
                                                                       4446, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 7981, 0, 3, 7621,
                                                                       4146, 7666, 2114, 2150,
                                                                       4506, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 8071, 0, 3, 7666,
                                                                       4176, 7711, 2150, 2186,
                                                                       4566, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 8161, 0, 3, 7711,
                                                                       4206, 7756, 2186, 2222,
                                                                       4626, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 8251, 0, 3, 7756,
                                                                       4236, 7801, 2222, 2258,
                                                                       4686, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 8341, 0, 3, 7801,
                                                                       4266, 7846, 2258, 2294,
                                                                       4746, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 8431, 0, 3, 7891,
                                                                       4446, 7981, 2366, 2426,
                                                                       5006, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 8581, 0, 3, 7981,
                                                                       4506, 8071, 2426, 2486,
                                                                       5106, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 8731, 0, 3, 8071,
                                                                       4566, 8161, 2486, 2546,
                                                                       5206, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 8881, 0, 3, 8161,
                                                                       4626, 8251, 2546, 2606,
                                                                       5306, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 9031, 0, 3, 8251,
                                                                       4686, 8341, 2606, 2666,
                                                                       5406, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 9181, 0, 3, 8431,
                                                                       5006, 8581, 2786, 2876,
                                                                       5806, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 9406, 0, 3, 8581,
                                                                       5106, 8731, 2876, 2966,
                                                                       5956, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 9631, 0, 3, 8731,
                                                                       5206, 8881, 2966, 3056,
                                                                       6106, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 9856, 0, 3, 8881,
                                                                       5306, 9031, 3056, 3146,
                                                                       6256, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 10081, 0, 3, 9181,
                                                                       5806, 9406, 3326, 3452,
                                                                       6826, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 10396, 0, 3, 9406,
                                                                       5956, 9631, 3452, 3578,
                                                                       7036, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 10711, 0, 3, 9631,
                                                                       6106, 9856, 3578, 3704,
                                                                       7246, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11026, 3, 3956,
                                                                       3966, 7456, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11047, 3, 3966,
                                                                       3976, 7471, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11068, 3, 3976,
                                                                       3986, 7486, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11089, 3, 3986,
                                                                       3996, 7501, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11110, 3, 3996,
                                                                       4006, 7516, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11131, 3, 4006,
                                                                       4016, 7531, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11152, 3, 4016,
                                                                       4026, 7546, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11173, 3, 4026,
                                                                       4036, 7561, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 11194, 0, 3,
                                                                       11026, 7456, 11047, 4056,
                                                                       4086, 7576, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 11257, 0, 3,
                                                                       11047, 7471, 11068, 4086,
                                                                       4116, 7621, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 11320, 0, 3,
                                                                       11068, 7486, 11089, 4116,
                                                                       4146, 7666, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 11383, 0, 3,
                                                                       11089, 7501, 11110, 4146,
                                                                       4176, 7711, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 11446, 0, 3,
                                                                       11110, 7516, 11131, 4176,
                                                                       4206, 7756, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 11509, 0, 3,
                                                                       11131, 7531, 11152, 4206,
                                                                       4236, 7801, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 11572, 0, 3,
                                                                       11152, 7546, 11173, 4236,
                                                                       4266, 7846, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 11635, 0, 3,
                                                                       11194, 7576, 11257, 4326,
                                                                       4386, 7891, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 11761, 0, 3,
                                                                       11257, 7621, 11320, 4386,
                                                                       4446, 7981, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 11887, 0, 3,
                                                                       11320, 7666, 11383, 4446,
                                                                       4506, 8071, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 12013, 0, 3,
                                                                       11383, 7711, 11446, 4506,
                                                                       4566, 8161, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 12139, 0, 3,
                                                                       11446, 7756, 11509, 4566,
                                                                       4626, 8251, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 12265, 0, 3,
                                                                       11509, 7801, 11572, 4626,
                                                                       4686, 8341, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 12391, 0, 3,
                                                                       11635, 7891, 11761, 4806,
                                                                       4906, 8431, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 12601, 0, 3,
                                                                       11761, 7981, 11887, 4906,
                                                                       5006, 8581, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 12811, 0, 3,
                                                                       11887, 8071, 12013, 5006,
                                                                       5106, 8731, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 13021, 0, 3,
                                                                       12013, 8161, 12139, 5106,
                                                                       5206, 8881, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 13231, 0, 3,
                                                                       12139, 8251, 12265, 5206,
                                                                       5306, 9031, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 13441, 0, 3,
                                                                       12391, 8431, 12601, 5506,
                                                                       5656, 9181, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 13756, 0, 3,
                                                                       12601, 8581, 12811, 5656,
                                                                       5806, 9406, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 14071, 0, 3,
                                                                       12811, 8731, 13021, 5806,
                                                                       5956, 9631, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 14386, 0, 3,
                                                                       13021, 8881, 13231, 5956,
                                                                       6106, 9856, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 14701, 0, 3,
                                                                       13441, 9181, 13756, 6406,
                                                                       6616, 10081, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 15142, 0, 3,
                                                                       13756, 9406, 14071, 6616,
                                                                       6826, 10396, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 15583, 0, 3,
                                                                       14071, 9631, 14386, 6826,
                                                                       7036, 10711, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16024, 3, 7456,
                                                                       7471, 11068, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16052, 3, 7471,
                                                                       7486, 11089, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16080, 3, 7486,
                                                                       7501, 11110, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16108, 3, 7501,
                                                                       7516, 11131, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16136, 3, 7516,
                                                                       7531, 11152, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16164, 3, 7531,
                                                                       7546, 11173, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 16192, 0, 3,
                                                                       16024, 11068, 16052, 7576,
                                                                       7621, 11320, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 16276, 0, 3,
                                                                       16052, 11089, 16080, 7621,
                                                                       7666, 11383, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 16360, 0, 3,
                                                                       16080, 11110, 16108, 7666,
                                                                       7711, 11446, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 16444, 0, 3,
                                                                       16108, 11131, 16136, 7711,
                                                                       7756, 11509, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 16528, 0, 3,
                                                                       16136, 11152, 16164, 7756,
                                                                       7801, 11572, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 16612, 0, 3,
                                                                       16192, 11320, 16276, 7891,
                                                                       7981, 11887, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 16780, 0, 3,
                                                                       16276, 11383, 16360, 7981,
                                                                       8071, 12013, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 16948, 0, 3,
                                                                       16360, 11446, 16444, 8071,
                                                                       8161, 12139, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 17116, 0, 3,
                                                                       16444, 11509, 16528, 8161,
                                                                       8251, 12265, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 17284, 0, 3,
                                                                       16612, 11887, 16780, 8431,
                                                                       8581, 12811, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 17564, 0, 3,
                                                                       16780, 12013, 16948, 8581,
                                                                       8731, 13021, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 17844, 0, 3,
                                                                       16948, 12139, 17116, 8731,
                                                                       8881, 13231, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 18124, 0, 3,
                                                                       17284, 12811, 17564, 9181,
                                                                       9406, 14071, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 18544, 0, 3,
                                                                       17564, 13021, 17844, 9406,
                                                                       9631, 14386, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 18964, 0, 3,
                                                                       18124, 14071, 18544,
                                                                       10081, 10396, 15583,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 19552, 3, 11026,
                                                                       11047, 16024, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 19588, 3, 11047,
                                                                       11068, 16052, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 19624, 3, 11068,
                                                                       11089, 16080, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 19660, 3, 11089,
                                                                       11110, 16108, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 19696, 3, 11110,
                                                                       11131, 16136, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 19732, 3, 11131,
                                                                       11152, 16164, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 19768, 0, 3,
                                                                       19552, 16024, 19588,
                                                                       11194, 11257, 16192,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 19876, 0, 3,
                                                                       19588, 16052, 19624,
                                                                       11257, 11320, 16276,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 19984, 0, 3,
                                                                       19624, 16080, 19660,
                                                                       11320, 11383, 16360,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 20092, 0, 3,
                                                                       19660, 16108, 19696,
                                                                       11383, 11446, 16444,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 20200, 0, 3,
                                                                       19696, 16136, 19732,
                                                                       11446, 11509, 16528,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 20308, 0, 3,
                                                                       19768, 16192, 19876,
                                                                       11635, 11761, 16612,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 20524, 0, 3,
                                                                       19876, 16276, 19984,
                                                                       11761, 11887, 16780,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 20740, 0, 3,
                                                                       19984, 16360, 20092,
                                                                       11887, 12013, 16948,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 20956, 0, 3,
                                                                       20092, 16444, 20200,
                                                                       12013, 12139, 17116,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 21172, 0, 3,
                                                                       20308, 16612, 20524,
                                                                       12391, 12601, 17284,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 21532, 0, 3,
                                                                       20524, 16780, 20740,
                                                                       12601, 12811, 17564,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 21892, 0, 3,
                                                                       20740, 16948, 20956,
                                                                       12811, 13021, 17844,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 22252, 0, 3,
                                                                       21172, 17284, 21532,
                                                                       13441, 13756, 18124,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 22792, 0, 3,
                                                                       21532, 17564, 21892,
                                                                       13756, 14071, 18544,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 23332, 0, 3,
                                                                       22252, 18124, 22792,
                                                                       14701, 15142, 18964,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 24088, 21172, 360, ncols);

                    simdfunc::contract_primitives(buffer, 24598, 22252, 540, ncols);

                    simdfunc::contract_primitives(buffer, 25363, 23332, 756, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 24448, 24088, 10, 1, nmax);

        simdtrf::transform_k_inner(buffer, 25138, 24598, 15, 1, nmax);

        simdtrf::transform_k_inner(buffer, 26119, 25363, 21, 1, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 26434, 24448, 25138, 15,
                                             nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 26884, 25138, 26119, 15,
                                             nmax);

        simdtrf::compute_hrr_fd_out_of_first(buffer, coordinates, 27559, 26434, 26884, 15,
                                             nmax);

        simdtrf::transform_d_inner(buffer, 28459, 27559, 10, 15, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 28459, 75, nmax);
    }

    for (size_t m = 0; m < 525; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
