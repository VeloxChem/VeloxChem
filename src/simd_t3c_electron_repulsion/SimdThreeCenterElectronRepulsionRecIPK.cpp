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


#include "SimdThreeCenterElectronRepulsionRecIPK.hpp"

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
#include "SimdTransferIP.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformK.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_ipk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_ipk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 64592, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 585 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 64592, 58808, 2724, dimensions);

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
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14},
                                                        ncols, fj, 6, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 22, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 25, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 28, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 31, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 34, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 37, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 40, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 43, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 46, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 49, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 52, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 55, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 58, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 61, 0, 3, 8, 9,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 67, 0, 3, 9, 10,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 73, 0, 3, 10, 11,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 79, 0, 3, 11, 12,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 85, 0, 3, 12, 13,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 91, 0, 3, 13, 14,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 97, 0, 3, 14, 15,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 103, 0, 3, 15, 16,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 109, 0, 3, 16, 17,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 115, 0, 3, 17, 18,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 121, 0, 3, 18, 19,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 127, 0, 3, 19, 20,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 133, 0, 3, 22, 25,
                                                                       61, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 143, 0, 3, 25, 28,
                                                                       67, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 153, 0, 3, 28, 31,
                                                                       73, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 163, 0, 3, 31, 34,
                                                                       79, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 173, 0, 3, 34, 37,
                                                                       85, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 183, 0, 3, 37, 40,
                                                                       91, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 193, 0, 3, 40, 43,
                                                                       97, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 203, 0, 3, 43, 46,
                                                                       103, 109, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 213, 0, 3, 46, 49,
                                                                       109, 115, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 223, 0, 3, 49, 52,
                                                                       115, 121, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 233, 0, 3, 52, 55,
                                                                       121, 127, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 243, 0, 3, 61, 67,
                                                                       133, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 258, 0, 3, 67, 73,
                                                                       143, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 273, 0, 3, 73, 79,
                                                                       153, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 288, 0, 3, 79, 85,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 303, 0, 3, 85, 91,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 318, 0, 3, 91, 97,
                                                                       183, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 333, 0, 3, 97,
                                                                       103, 193, 203, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 348, 0, 3, 103,
                                                                       109, 203, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 363, 0, 3, 109,
                                                                       115, 213, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 378, 0, 3, 115,
                                                                       121, 223, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 393, 0, 3, 133,
                                                                       143, 243, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 414, 0, 3, 143,
                                                                       153, 258, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 435, 0, 3, 153,
                                                                       163, 273, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 456, 0, 3, 163,
                                                                       173, 288, 303, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 477, 0, 3, 173,
                                                                       183, 303, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 498, 0, 3, 183,
                                                                       193, 318, 333, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 519, 0, 3, 193,
                                                                       203, 333, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 540, 0, 3, 203,
                                                                       213, 348, 363, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 561, 0, 3, 213,
                                                                       223, 363, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 582, 0, 3, 243,
                                                                       258, 393, 414, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 610, 0, 3, 258,
                                                                       273, 414, 435, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 638, 0, 3, 273,
                                                                       288, 435, 456, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 666, 0, 3, 288,
                                                                       303, 456, 477, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 694, 0, 3, 303,
                                                                       318, 477, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 722, 0, 3, 318,
                                                                       333, 498, 519, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 750, 0, 3, 333,
                                                                       348, 519, 540, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 778, 0, 3, 348,
                                                                       363, 540, 561, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 806, 0, 3, 393,
                                                                       414, 582, 610, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 842, 0, 3, 414,
                                                                       435, 610, 638, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 878, 0, 3, 435,
                                                                       456, 638, 666, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 914, 0, 3, 456,
                                                                       477, 666, 694, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 950, 0, 3, 477,
                                                                       498, 694, 722, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 986, 0, 3, 498,
                                                                       519, 722, 750, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1022, 0, 3, 519,
                                                                       540, 750, 778, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1058, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1061, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1064, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1067, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1070, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1073, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1076, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1079, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1082, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1085, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1088, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1091, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1094, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1097, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1100, 3, 8, 22,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1109, 3, 9, 25,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1118, 3, 10, 28,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1127, 3, 11, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1136, 3, 12, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1145, 3, 13, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1154, 3, 14, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1163, 3, 15, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1172, 3, 16, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1181, 3, 17, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1190, 3, 18, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1199, 3, 19, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1208, 3, 20, 58,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1217, 3, 22, 61,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1235, 3, 25, 67,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1253, 3, 28, 73,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1271, 3, 31, 79,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1289, 3, 34, 85,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1307, 3, 37, 91,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1325, 3, 40, 97,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1343, 3, 43, 103,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1361, 3, 46, 109,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1379, 3, 49, 115,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1397, 3, 52, 121,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1415, 3, 55, 127,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1433, 3, 61, 133,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1463, 3, 67, 143,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1493, 3, 73, 153,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1523, 3, 79, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1553, 3, 85, 173,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1583, 3, 91, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1613, 3, 97, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1643, 3, 103, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1673, 3, 109, 213,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1703, 3, 115, 223,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1733, 3, 121, 233,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1763, 3, 133, 243,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1808, 3, 143, 258,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1853, 3, 153, 273,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1898, 3, 163, 288,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1943, 3, 173, 303,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1988, 3, 183, 318,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2033, 3, 193, 333,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2078, 3, 203, 348,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2123, 3, 213, 363,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2168, 3, 223, 378,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2213, 3, 243, 393,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2276, 3, 258, 414,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2339, 3, 273, 435,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2402, 3, 288, 456,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2465, 3, 303, 477,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2528, 3, 318, 498,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2591, 3, 333, 519,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2654, 3, 348, 540,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2717, 3, 363, 561,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2780, 3, 393, 582,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2864, 3, 414, 610,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2948, 3, 435, 638,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3032, 3, 456, 666,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3116, 3, 477, 694,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3200, 3, 498, 722,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3284, 3, 519, 750,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3368, 3, 540, 778,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3452, 3, 582, 806,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3560, 3, 610, 842,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3668, 3, 638, 878,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3776, 3, 666, 914,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3884, 3, 694, 950,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3992, 3, 722, 986,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4100, 3, 750,
                                                                       1022, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4208, 3, 8, 9,
                                                                       1064, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4214, 3, 9, 10,
                                                                       1067, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4220, 3, 10, 11,
                                                                       1070, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4226, 3, 11, 12,
                                                                       1073, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4232, 3, 12, 13,
                                                                       1076, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4238, 3, 13, 14,
                                                                       1079, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4244, 3, 14, 15,
                                                                       1082, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4250, 3, 15, 16,
                                                                       1085, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4256, 3, 16, 17,
                                                                       1088, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4262, 3, 17, 18,
                                                                       1091, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4268, 3, 18, 19,
                                                                       1094, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4274, 3, 19, 20,
                                                                       1097, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4280, 0, 3, 4208,
                                                                       1064, 4214, 1118, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4298, 0, 3, 4214,
                                                                       1067, 4220, 1127, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4316, 0, 3, 4220,
                                                                       1070, 4226, 1136, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4334, 0, 3, 4226,
                                                                       1073, 4232, 1145, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4352, 0, 3, 4232,
                                                                       1076, 4238, 1154, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4370, 0, 3, 4238,
                                                                       1079, 4244, 1163, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4388, 0, 3, 4244,
                                                                       1082, 4250, 1172, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4406, 0, 3, 4250,
                                                                       1085, 4256, 1181, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4424, 0, 3, 4256,
                                                                       1088, 4262, 1190, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4442, 0, 3, 4262,
                                                                       1091, 4268, 1199, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4460, 0, 3, 4268,
                                                                       1094, 4274, 1208, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4478, 0, 3, 4280,
                                                                       1118, 4298, 61, 67, 1253,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4514, 0, 3, 4298,
                                                                       1127, 4316, 67, 73, 1271,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4550, 0, 3, 4316,
                                                                       1136, 4334, 73, 79, 1289,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4586, 0, 3, 4334,
                                                                       1145, 4352, 79, 85, 1307,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4622, 0, 3, 4352,
                                                                       1154, 4370, 85, 91, 1325,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4658, 0, 3, 4370,
                                                                       1163, 4388, 91, 97, 1343,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4694, 0, 3, 4388,
                                                                       1172, 4406, 97, 103, 1361,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4730, 0, 3, 4406,
                                                                       1181, 4424, 103, 109,
                                                                       1379, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4766, 0, 3, 4424,
                                                                       1190, 4442, 109, 115,
                                                                       1397, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4802, 0, 3, 4442,
                                                                       1199, 4460, 115, 121,
                                                                       1415, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4838, 0, 3, 4478,
                                                                       1253, 4514, 133, 143,
                                                                       1493, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4898, 0, 3, 4514,
                                                                       1271, 4550, 143, 153,
                                                                       1523, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4958, 0, 3, 4550,
                                                                       1289, 4586, 153, 163,
                                                                       1553, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5018, 0, 3, 4586,
                                                                       1307, 4622, 163, 173,
                                                                       1583, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5078, 0, 3, 4622,
                                                                       1325, 4658, 173, 183,
                                                                       1613, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5138, 0, 3, 4658,
                                                                       1343, 4694, 183, 193,
                                                                       1643, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5198, 0, 3, 4694,
                                                                       1361, 4730, 193, 203,
                                                                       1673, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5258, 0, 3, 4730,
                                                                       1379, 4766, 203, 213,
                                                                       1703, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5318, 0, 3, 4766,
                                                                       1397, 4802, 213, 223,
                                                                       1733, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5378, 0, 3, 4838,
                                                                       1493, 4898, 243, 258,
                                                                       1853, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5468, 0, 3, 4898,
                                                                       1523, 4958, 258, 273,
                                                                       1898, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5558, 0, 3, 4958,
                                                                       1553, 5018, 273, 288,
                                                                       1943, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5648, 0, 3, 5018,
                                                                       1583, 5078, 288, 303,
                                                                       1988, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5738, 0, 3, 5078,
                                                                       1613, 5138, 303, 318,
                                                                       2033, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5828, 0, 3, 5138,
                                                                       1643, 5198, 318, 333,
                                                                       2078, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5918, 0, 3, 5198,
                                                                       1673, 5258, 333, 348,
                                                                       2123, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6008, 0, 3, 5258,
                                                                       1703, 5318, 348, 363,
                                                                       2168, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6098, 0, 3, 5378,
                                                                       1853, 5468, 393, 414,
                                                                       2339, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6224, 0, 3, 5468,
                                                                       1898, 5558, 414, 435,
                                                                       2402, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6350, 0, 3, 5558,
                                                                       1943, 5648, 435, 456,
                                                                       2465, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6476, 0, 3, 5648,
                                                                       1988, 5738, 456, 477,
                                                                       2528, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6602, 0, 3, 5738,
                                                                       2033, 5828, 477, 498,
                                                                       2591, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6728, 0, 3, 5828,
                                                                       2078, 5918, 498, 519,
                                                                       2654, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6854, 0, 3, 5918,
                                                                       2123, 6008, 519, 540,
                                                                       2717, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 6980, 0, 3, 6098,
                                                                       2339, 6224, 582, 610,
                                                                       2948, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7148, 0, 3, 6224,
                                                                       2402, 6350, 610, 638,
                                                                       3032, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7316, 0, 3, 6350,
                                                                       2465, 6476, 638, 666,
                                                                       3116, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7484, 0, 3, 6476,
                                                                       2528, 6602, 666, 694,
                                                                       3200, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7652, 0, 3, 6602,
                                                                       2591, 6728, 694, 722,
                                                                       3284, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7820, 0, 3, 6728,
                                                                       2654, 6854, 722, 750,
                                                                       3368, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 7988, 0, 3, 6980,
                                                                       2948, 7148, 806, 842,
                                                                       3668, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 8204, 0, 3, 7148,
                                                                       3032, 7316, 842, 878,
                                                                       3776, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 8420, 0, 3, 7316,
                                                                       3116, 7484, 878, 914,
                                                                       3884, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 8636, 0, 3, 7484,
                                                                       3200, 7652, 914, 950,
                                                                       3992, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 8852, 0, 3, 7652,
                                                                       3284, 7820, 950, 986,
                                                                       4100, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9068, 3, 1058,
                                                                       1061, 4208, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9078, 3, 1061,
                                                                       1064, 4214, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9088, 3, 1064,
                                                                       1067, 4220, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9098, 3, 1067,
                                                                       1070, 4226, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9108, 3, 1070,
                                                                       1073, 4232, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9118, 3, 1073,
                                                                       1076, 4238, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9128, 3, 1076,
                                                                       1079, 4244, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9138, 3, 1079,
                                                                       1082, 4250, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9148, 3, 1082,
                                                                       1085, 4256, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9158, 3, 1085,
                                                                       1088, 4262, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9168, 3, 1088,
                                                                       1091, 4268, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9178, 3, 1091,
                                                                       1094, 4274, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9188, 0, 3, 9068,
                                                                       4208, 9078, 1100, 1109,
                                                                       4280, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9218, 0, 3, 9078,
                                                                       4214, 9088, 1109, 1118,
                                                                       4298, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9248, 0, 3, 9088,
                                                                       4220, 9098, 1118, 1127,
                                                                       4316, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9278, 0, 3, 9098,
                                                                       4226, 9108, 1127, 1136,
                                                                       4334, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9308, 0, 3, 9108,
                                                                       4232, 9118, 1136, 1145,
                                                                       4352, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9338, 0, 3, 9118,
                                                                       4238, 9128, 1145, 1154,
                                                                       4370, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9368, 0, 3, 9128,
                                                                       4244, 9138, 1154, 1163,
                                                                       4388, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9398, 0, 3, 9138,
                                                                       4250, 9148, 1163, 1172,
                                                                       4406, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9428, 0, 3, 9148,
                                                                       4256, 9158, 1172, 1181,
                                                                       4424, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9458, 0, 3, 9158,
                                                                       4262, 9168, 1181, 1190,
                                                                       4442, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9488, 0, 3, 9168,
                                                                       4268, 9178, 1190, 1199,
                                                                       4460, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9518, 0, 3, 9188,
                                                                       4280, 9218, 1217, 1235,
                                                                       4478, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9578, 0, 3, 9218,
                                                                       4298, 9248, 1235, 1253,
                                                                       4514, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9638, 0, 3, 9248,
                                                                       4316, 9278, 1253, 1271,
                                                                       4550, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9698, 0, 3, 9278,
                                                                       4334, 9308, 1271, 1289,
                                                                       4586, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9758, 0, 3, 9308,
                                                                       4352, 9338, 1289, 1307,
                                                                       4622, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9818, 0, 3, 9338,
                                                                       4370, 9368, 1307, 1325,
                                                                       4658, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9878, 0, 3, 9368,
                                                                       4388, 9398, 1325, 1343,
                                                                       4694, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9938, 0, 3, 9398,
                                                                       4406, 9428, 1343, 1361,
                                                                       4730, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9998, 0, 3, 9428,
                                                                       4424, 9458, 1361, 1379,
                                                                       4766, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10058, 0, 3, 9458,
                                                                       4442, 9488, 1379, 1397,
                                                                       4802, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10118, 0, 3, 9518,
                                                                       4478, 9578, 1433, 1463,
                                                                       4838, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10218, 0, 3, 9578,
                                                                       4514, 9638, 1463, 1493,
                                                                       4898, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10318, 0, 3, 9638,
                                                                       4550, 9698, 1493, 1523,
                                                                       4958, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10418, 0, 3, 9698,
                                                                       4586, 9758, 1523, 1553,
                                                                       5018, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10518, 0, 3, 9758,
                                                                       4622, 9818, 1553, 1583,
                                                                       5078, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10618, 0, 3, 9818,
                                                                       4658, 9878, 1583, 1613,
                                                                       5138, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10718, 0, 3, 9878,
                                                                       4694, 9938, 1613, 1643,
                                                                       5198, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10818, 0, 3, 9938,
                                                                       4730, 9998, 1643, 1673,
                                                                       5258, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10918, 0, 3, 9998,
                                                                       4766, 10058, 1673, 1703,
                                                                       5318, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11018, 0, 3,
                                                                       10118, 4838, 10218, 1763,
                                                                       1808, 5378, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11168, 0, 3,
                                                                       10218, 4898, 10318, 1808,
                                                                       1853, 5468, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11318, 0, 3,
                                                                       10318, 4958, 10418, 1853,
                                                                       1898, 5558, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11468, 0, 3,
                                                                       10418, 5018, 10518, 1898,
                                                                       1943, 5648, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11618, 0, 3,
                                                                       10518, 5078, 10618, 1943,
                                                                       1988, 5738, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11768, 0, 3,
                                                                       10618, 5138, 10718, 1988,
                                                                       2033, 5828, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11918, 0, 3,
                                                                       10718, 5198, 10818, 2033,
                                                                       2078, 5918, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 12068, 0, 3,
                                                                       10818, 5258, 10918, 2078,
                                                                       2123, 6008, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 12218, 0, 3,
                                                                       11018, 5378, 11168, 2213,
                                                                       2276, 6098, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 12428, 0, 3,
                                                                       11168, 5468, 11318, 2276,
                                                                       2339, 6224, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 12638, 0, 3,
                                                                       11318, 5558, 11468, 2339,
                                                                       2402, 6350, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 12848, 0, 3,
                                                                       11468, 5648, 11618, 2402,
                                                                       2465, 6476, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 13058, 0, 3,
                                                                       11618, 5738, 11768, 2465,
                                                                       2528, 6602, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 13268, 0, 3,
                                                                       11768, 5828, 11918, 2528,
                                                                       2591, 6728, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 13478, 0, 3,
                                                                       11918, 5918, 12068, 2591,
                                                                       2654, 6854, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 13688, 0, 3,
                                                                       12218, 6098, 12428, 2780,
                                                                       2864, 6980, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 13968, 0, 3,
                                                                       12428, 6224, 12638, 2864,
                                                                       2948, 7148, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 14248, 0, 3,
                                                                       12638, 6350, 12848, 2948,
                                                                       3032, 7316, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 14528, 0, 3,
                                                                       12848, 6476, 13058, 3032,
                                                                       3116, 7484, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 14808, 0, 3,
                                                                       13058, 6602, 13268, 3116,
                                                                       3200, 7652, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 15088, 0, 3,
                                                                       13268, 6728, 13478, 3200,
                                                                       3284, 7820, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 15368, 0, 3,
                                                                       13688, 6980, 13968, 3452,
                                                                       3560, 7988, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 15728, 0, 3,
                                                                       13968, 7148, 14248, 3560,
                                                                       3668, 8204, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 16088, 0, 3,
                                                                       14248, 7316, 14528, 3668,
                                                                       3776, 8420, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 16448, 0, 3,
                                                                       14528, 7484, 14808, 3776,
                                                                       3884, 8636, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 16808, 0, 3,
                                                                       14808, 7652, 15088, 3884,
                                                                       3992, 8852, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17168, 3, 4208,
                                                                       4214, 9088, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17183, 3, 4214,
                                                                       4220, 9098, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17198, 3, 4220,
                                                                       4226, 9108, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17213, 3, 4226,
                                                                       4232, 9118, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17228, 3, 4232,
                                                                       4238, 9128, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17243, 3, 4238,
                                                                       4244, 9138, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17258, 3, 4244,
                                                                       4250, 9148, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17273, 3, 4250,
                                                                       4256, 9158, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17288, 3, 4256,
                                                                       4262, 9168, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17303, 3, 4262,
                                                                       4268, 9178, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17318, 0, 3,
                                                                       17168, 9088, 17183, 4280,
                                                                       4298, 9248, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17363, 0, 3,
                                                                       17183, 9098, 17198, 4298,
                                                                       4316, 9278, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17408, 0, 3,
                                                                       17198, 9108, 17213, 4316,
                                                                       4334, 9308, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17453, 0, 3,
                                                                       17213, 9118, 17228, 4334,
                                                                       4352, 9338, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17498, 0, 3,
                                                                       17228, 9128, 17243, 4352,
                                                                       4370, 9368, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17543, 0, 3,
                                                                       17243, 9138, 17258, 4370,
                                                                       4388, 9398, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17588, 0, 3,
                                                                       17258, 9148, 17273, 4388,
                                                                       4406, 9428, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17633, 0, 3,
                                                                       17273, 9158, 17288, 4406,
                                                                       4424, 9458, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17678, 0, 3,
                                                                       17288, 9168, 17303, 4424,
                                                                       4442, 9488, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 17723, 0, 3,
                                                                       17318, 9248, 17363, 4478,
                                                                       4514, 9638, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 17813, 0, 3,
                                                                       17363, 9278, 17408, 4514,
                                                                       4550, 9698, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 17903, 0, 3,
                                                                       17408, 9308, 17453, 4550,
                                                                       4586, 9758, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 17993, 0, 3,
                                                                       17453, 9338, 17498, 4586,
                                                                       4622, 9818, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 18083, 0, 3,
                                                                       17498, 9368, 17543, 4622,
                                                                       4658, 9878, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 18173, 0, 3,
                                                                       17543, 9398, 17588, 4658,
                                                                       4694, 9938, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 18263, 0, 3,
                                                                       17588, 9428, 17633, 4694,
                                                                       4730, 9998, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 18353, 0, 3,
                                                                       17633, 9458, 17678, 4730,
                                                                       4766, 10058, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 18443, 0, 3,
                                                                       17723, 9638, 17813, 4838,
                                                                       4898, 10318, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 18593, 0, 3,
                                                                       17813, 9698, 17903, 4898,
                                                                       4958, 10418, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 18743, 0, 3,
                                                                       17903, 9758, 17993, 4958,
                                                                       5018, 10518, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 18893, 0, 3,
                                                                       17993, 9818, 18083, 5018,
                                                                       5078, 10618, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 19043, 0, 3,
                                                                       18083, 9878, 18173, 5078,
                                                                       5138, 10718, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 19193, 0, 3,
                                                                       18173, 9938, 18263, 5138,
                                                                       5198, 10818, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 19343, 0, 3,
                                                                       18263, 9998, 18353, 5198,
                                                                       5258, 10918, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 19493, 0, 3,
                                                                       18443, 10318, 18593, 5378,
                                                                       5468, 11318, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 19718, 0, 3,
                                                                       18593, 10418, 18743, 5468,
                                                                       5558, 11468, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 19943, 0, 3,
                                                                       18743, 10518, 18893, 5558,
                                                                       5648, 11618, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 20168, 0, 3,
                                                                       18893, 10618, 19043, 5648,
                                                                       5738, 11768, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 20393, 0, 3,
                                                                       19043, 10718, 19193, 5738,
                                                                       5828, 11918, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 20618, 0, 3,
                                                                       19193, 10818, 19343, 5828,
                                                                       5918, 12068, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 20843, 0, 3,
                                                                       19493, 11318, 19718, 6098,
                                                                       6224, 12638, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 21158, 0, 3,
                                                                       19718, 11468, 19943, 6224,
                                                                       6350, 12848, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 21473, 0, 3,
                                                                       19943, 11618, 20168, 6350,
                                                                       6476, 13058, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 21788, 0, 3,
                                                                       20168, 11768, 20393, 6476,
                                                                       6602, 13268, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 22103, 0, 3,
                                                                       20393, 11918, 20618, 6602,
                                                                       6728, 13478, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 22418, 0, 3,
                                                                       20843, 12638, 21158, 6980,
                                                                       7148, 14248, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 22838, 0, 3,
                                                                       21158, 12848, 21473, 7148,
                                                                       7316, 14528, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 23258, 0, 3,
                                                                       21473, 13058, 21788, 7316,
                                                                       7484, 14808, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 23678, 0, 3,
                                                                       21788, 13268, 22103, 7484,
                                                                       7652, 15088, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 24098, 0, 3,
                                                                       22418, 14248, 22838, 7988,
                                                                       8204, 16088, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 24638, 0, 3,
                                                                       22838, 14528, 23258, 8204,
                                                                       8420, 16448, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 25178, 0, 3,
                                                                       23258, 14808, 23678, 8420,
                                                                       8636, 16808, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25718, 3, 9068,
                                                                       9078, 17168, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25739, 3, 9078,
                                                                       9088, 17183, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25760, 3, 9088,
                                                                       9098, 17198, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25781, 3, 9098,
                                                                       9108, 17213, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25802, 3, 9108,
                                                                       9118, 17228, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25823, 3, 9118,
                                                                       9128, 17243, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25844, 3, 9128,
                                                                       9138, 17258, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25865, 3, 9138,
                                                                       9148, 17273, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25886, 3, 9148,
                                                                       9158, 17288, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25907, 3, 9158,
                                                                       9168, 17303, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 25928, 0, 3,
                                                                       25718, 17168, 25739, 9188,
                                                                       9218, 17318, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 25991, 0, 3,
                                                                       25739, 17183, 25760, 9218,
                                                                       9248, 17363, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 26054, 0, 3,
                                                                       25760, 17198, 25781, 9248,
                                                                       9278, 17408, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 26117, 0, 3,
                                                                       25781, 17213, 25802, 9278,
                                                                       9308, 17453, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 26180, 0, 3,
                                                                       25802, 17228, 25823, 9308,
                                                                       9338, 17498, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 26243, 0, 3,
                                                                       25823, 17243, 25844, 9338,
                                                                       9368, 17543, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 26306, 0, 3,
                                                                       25844, 17258, 25865, 9368,
                                                                       9398, 17588, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 26369, 0, 3,
                                                                       25865, 17273, 25886, 9398,
                                                                       9428, 17633, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 26432, 0, 3,
                                                                       25886, 17288, 25907, 9428,
                                                                       9458, 17678, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 26495, 0, 3,
                                                                       25928, 17318, 25991, 9518,
                                                                       9578, 17723, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 26621, 0, 3,
                                                                       25991, 17363, 26054, 9578,
                                                                       9638, 17813, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 26747, 0, 3,
                                                                       26054, 17408, 26117, 9638,
                                                                       9698, 17903, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 26873, 0, 3,
                                                                       26117, 17453, 26180, 9698,
                                                                       9758, 17993, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 26999, 0, 3,
                                                                       26180, 17498, 26243, 9758,
                                                                       9818, 18083, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 27125, 0, 3,
                                                                       26243, 17543, 26306, 9818,
                                                                       9878, 18173, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 27251, 0, 3,
                                                                       26306, 17588, 26369, 9878,
                                                                       9938, 18263, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 27377, 0, 3,
                                                                       26369, 17633, 26432, 9938,
                                                                       9998, 18353, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 27503, 0, 3,
                                                                       26495, 17723, 26621,
                                                                       10118, 10218, 18443,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 27713, 0, 3,
                                                                       26621, 17813, 26747,
                                                                       10218, 10318, 18593,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 27923, 0, 3,
                                                                       26747, 17903, 26873,
                                                                       10318, 10418, 18743,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 28133, 0, 3,
                                                                       26873, 17993, 26999,
                                                                       10418, 10518, 18893,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 28343, 0, 3,
                                                                       26999, 18083, 27125,
                                                                       10518, 10618, 19043,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 28553, 0, 3,
                                                                       27125, 18173, 27251,
                                                                       10618, 10718, 19193,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 28763, 0, 3,
                                                                       27251, 18263, 27377,
                                                                       10718, 10818, 19343,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 28973, 0, 3,
                                                                       27503, 18443, 27713,
                                                                       11018, 11168, 19493,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 29288, 0, 3,
                                                                       27713, 18593, 27923,
                                                                       11168, 11318, 19718,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 29603, 0, 3,
                                                                       27923, 18743, 28133,
                                                                       11318, 11468, 19943,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 29918, 0, 3,
                                                                       28133, 18893, 28343,
                                                                       11468, 11618, 20168,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 30233, 0, 3,
                                                                       28343, 19043, 28553,
                                                                       11618, 11768, 20393,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 30548, 0, 3,
                                                                       28553, 19193, 28763,
                                                                       11768, 11918, 20618,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 30863, 0, 3,
                                                                       28973, 19493, 29288,
                                                                       12218, 12428, 20843,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 31304, 0, 3,
                                                                       29288, 19718, 29603,
                                                                       12428, 12638, 21158,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 31745, 0, 3,
                                                                       29603, 19943, 29918,
                                                                       12638, 12848, 21473,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 32186, 0, 3,
                                                                       29918, 20168, 30233,
                                                                       12848, 13058, 21788,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 32627, 0, 3,
                                                                       30233, 20393, 30548,
                                                                       13058, 13268, 22103,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 33068, 0, 3,
                                                                       30863, 20843, 31304,
                                                                       13688, 13968, 22418,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 33656, 0, 3,
                                                                       31304, 21158, 31745,
                                                                       13968, 14248, 22838,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 34244, 0, 3,
                                                                       31745, 21473, 32186,
                                                                       14248, 14528, 23258,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 34832, 0, 3,
                                                                       32186, 21788, 32627,
                                                                       14528, 14808, 23678,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 35420, 0, 3,
                                                                       33068, 22418, 33656,
                                                                       15368, 15728, 24098,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 36176, 0, 3,
                                                                       33656, 22838, 34244,
                                                                       15728, 16088, 24638,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 36932, 0, 3,
                                                                       34244, 23258, 34832,
                                                                       16088, 16448, 25178,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37688, 3, 17168,
                                                                       17183, 25760, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37716, 3, 17183,
                                                                       17198, 25781, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37744, 3, 17198,
                                                                       17213, 25802, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37772, 3, 17213,
                                                                       17228, 25823, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37800, 3, 17228,
                                                                       17243, 25844, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37828, 3, 17243,
                                                                       17258, 25865, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37856, 3, 17258,
                                                                       17273, 25886, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37884, 3, 17273,
                                                                       17288, 25907, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 37912, 0, 3,
                                                                       37688, 25760, 37716,
                                                                       17318, 17363, 26054,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 37996, 0, 3,
                                                                       37716, 25781, 37744,
                                                                       17363, 17408, 26117,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 38080, 0, 3,
                                                                       37744, 25802, 37772,
                                                                       17408, 17453, 26180,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 38164, 0, 3,
                                                                       37772, 25823, 37800,
                                                                       17453, 17498, 26243,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 38248, 0, 3,
                                                                       37800, 25844, 37828,
                                                                       17498, 17543, 26306,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 38332, 0, 3,
                                                                       37828, 25865, 37856,
                                                                       17543, 17588, 26369,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 38416, 0, 3,
                                                                       37856, 25886, 37884,
                                                                       17588, 17633, 26432,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 38500, 0, 3,
                                                                       37912, 26054, 37996,
                                                                       17723, 17813, 26747,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 38668, 0, 3,
                                                                       37996, 26117, 38080,
                                                                       17813, 17903, 26873,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 38836, 0, 3,
                                                                       38080, 26180, 38164,
                                                                       17903, 17993, 26999,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 39004, 0, 3,
                                                                       38164, 26243, 38248,
                                                                       17993, 18083, 27125,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 39172, 0, 3,
                                                                       38248, 26306, 38332,
                                                                       18083, 18173, 27251,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 39340, 0, 3,
                                                                       38332, 26369, 38416,
                                                                       18173, 18263, 27377,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 39508, 0, 3,
                                                                       38500, 26747, 38668,
                                                                       18443, 18593, 27923,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 39788, 0, 3,
                                                                       38668, 26873, 38836,
                                                                       18593, 18743, 28133,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 40068, 0, 3,
                                                                       38836, 26999, 39004,
                                                                       18743, 18893, 28343,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 40348, 0, 3,
                                                                       39004, 27125, 39172,
                                                                       18893, 19043, 28553,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 40628, 0, 3,
                                                                       39172, 27251, 39340,
                                                                       19043, 19193, 28763,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 40908, 0, 3,
                                                                       39508, 27923, 39788,
                                                                       19493, 19718, 29603,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 41328, 0, 3,
                                                                       39788, 28133, 40068,
                                                                       19718, 19943, 29918,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 41748, 0, 3,
                                                                       40068, 28343, 40348,
                                                                       19943, 20168, 30233,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 42168, 0, 3,
                                                                       40348, 28553, 40628,
                                                                       20168, 20393, 30548,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 42588, 0, 3,
                                                                       40908, 29603, 41328,
                                                                       20843, 21158, 31745,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 43176, 0, 3,
                                                                       41328, 29918, 41748,
                                                                       21158, 21473, 32186,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 43764, 0, 3,
                                                                       41748, 30233, 42168,
                                                                       21473, 21788, 32627,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 44352, 0, 3,
                                                                       42588, 31745, 43176,
                                                                       22418, 22838, 34244,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 45136, 0, 3,
                                                                       43176, 32186, 43764,
                                                                       22838, 23258, 34832,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 45920, 0, 3,
                                                                       44352, 34244, 45136,
                                                                       24098, 24638, 36932,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 46928, 3, 25718,
                                                                       25739, 37688, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 46964, 3, 25739,
                                                                       25760, 37716, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 47000, 3, 25760,
                                                                       25781, 37744, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 47036, 3, 25781,
                                                                       25802, 37772, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 47072, 3, 25802,
                                                                       25823, 37800, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 47108, 3, 25823,
                                                                       25844, 37828, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 47144, 3, 25844,
                                                                       25865, 37856, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 47180, 3, 25865,
                                                                       25886, 37884, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 47216, 0, 3,
                                                                       46928, 37688, 46964,
                                                                       25928, 25991, 37912,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 47324, 0, 3,
                                                                       46964, 37716, 47000,
                                                                       25991, 26054, 37996,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 47432, 0, 3,
                                                                       47000, 37744, 47036,
                                                                       26054, 26117, 38080,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 47540, 0, 3,
                                                                       47036, 37772, 47072,
                                                                       26117, 26180, 38164,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 47648, 0, 3,
                                                                       47072, 37800, 47108,
                                                                       26180, 26243, 38248,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 47756, 0, 3,
                                                                       47108, 37828, 47144,
                                                                       26243, 26306, 38332,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 47864, 0, 3,
                                                                       47144, 37856, 47180,
                                                                       26306, 26369, 38416,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 47972, 0, 3,
                                                                       47216, 37912, 47324,
                                                                       26495, 26621, 38500,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 48188, 0, 3,
                                                                       47324, 37996, 47432,
                                                                       26621, 26747, 38668,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 48404, 0, 3,
                                                                       47432, 38080, 47540,
                                                                       26747, 26873, 38836,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 48620, 0, 3,
                                                                       47540, 38164, 47648,
                                                                       26873, 26999, 39004,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 48836, 0, 3,
                                                                       47648, 38248, 47756,
                                                                       26999, 27125, 39172,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 49052, 0, 3,
                                                                       47756, 38332, 47864,
                                                                       27125, 27251, 39340,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 49268, 0, 3,
                                                                       47972, 38500, 48188,
                                                                       27503, 27713, 39508,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 49628, 0, 3,
                                                                       48188, 38668, 48404,
                                                                       27713, 27923, 39788,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 49988, 0, 3,
                                                                       48404, 38836, 48620,
                                                                       27923, 28133, 40068,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 50348, 0, 3,
                                                                       48620, 39004, 48836,
                                                                       28133, 28343, 40348,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 50708, 0, 3,
                                                                       48836, 39172, 49052,
                                                                       28343, 28553, 40628,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 51068, 0, 3,
                                                                       49268, 39508, 49628,
                                                                       28973, 29288, 40908,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 51608, 0, 3,
                                                                       49628, 39788, 49988,
                                                                       29288, 29603, 41328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 52148, 0, 3,
                                                                       49988, 40068, 50348,
                                                                       29603, 29918, 41748,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 52688, 0, 3,
                                                                       50348, 40348, 50708,
                                                                       29918, 30233, 42168,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 53228, 0, 3,
                                                                       51068, 40908, 51608,
                                                                       30863, 31304, 42588,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 53984, 0, 3,
                                                                       51608, 41328, 52148,
                                                                       31304, 31745, 43176,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 54740, 0, 3,
                                                                       52148, 41748, 52688,
                                                                       31745, 32186, 43764,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 55496, 0, 3,
                                                                       53228, 42588, 53984,
                                                                       33068, 33656, 44352,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 56504, 0, 3,
                                                                       53984, 43176, 54740,
                                                                       33656, 34244, 45136,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 57512, 0, 3,
                                                                       55496, 44352, 56504,
                                                                       35420, 36176, 45920,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 58808, 55496, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 60236, 57512, 1296, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 59816, 58808, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 61532, 60236, 36, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 62072, 59816, 61532, 15,
                                             nmax);

        simdtrf::transform_p_inner(buffer, 63332, 62072, 28, 15, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 63332, 45, nmax);
    }

    for (size_t m = 0; m < 585; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
