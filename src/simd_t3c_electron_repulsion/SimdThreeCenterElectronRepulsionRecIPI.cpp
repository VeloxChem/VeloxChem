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


#include "SimdThreeCenterElectronRepulsionRecIPI.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_ipi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_ipi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 43576, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 507 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 43576, 38768, 2156, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 13,
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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1058, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1061, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1064, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1067, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1070, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1073, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1076, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1079, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1082, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1085, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1088, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1091, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1094, 3, 10, 28,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1103, 3, 11, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1112, 3, 12, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1121, 3, 13, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1130, 3, 14, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1139, 3, 15, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1148, 3, 16, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1157, 3, 17, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1166, 3, 18, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1175, 3, 19, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1184, 3, 20, 58,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1193, 3, 28, 73,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1211, 3, 31, 79,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1229, 3, 34, 85,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1247, 3, 37, 91,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1265, 3, 40, 97,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1283, 3, 43, 103,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1301, 3, 46, 109,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1319, 3, 49, 115,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1337, 3, 52, 121,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1355, 3, 55, 127,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1373, 3, 73, 153,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1403, 3, 79, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1433, 3, 85, 173,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1463, 3, 91, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1493, 3, 97, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1523, 3, 103, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1553, 3, 109, 213,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1583, 3, 115, 223,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1613, 3, 121, 233,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1643, 3, 153, 273,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1688, 3, 163, 288,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1733, 3, 173, 303,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1778, 3, 183, 318,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1823, 3, 193, 333,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1868, 3, 203, 348,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1913, 3, 213, 363,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1958, 3, 223, 378,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2003, 3, 273, 435,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2066, 3, 288, 456,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2129, 3, 303, 477,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2192, 3, 318, 498,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2255, 3, 333, 519,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2318, 3, 348, 540,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2381, 3, 363, 561,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2444, 3, 435, 638,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2528, 3, 456, 666,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2612, 3, 477, 694,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2696, 3, 498, 722,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2780, 3, 519, 750,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2864, 3, 540, 778,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 2948, 3, 638, 878,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3056, 3, 666, 914,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3164, 3, 694, 950,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3272, 3, 722, 986,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3380, 3, 750,
                                                                       1022, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3488, 3, 8, 9,
                                                                       1058, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3494, 3, 9, 10,
                                                                       1061, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3500, 3, 10, 11,
                                                                       1064, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3506, 3, 11, 12,
                                                                       1067, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3512, 3, 12, 13,
                                                                       1070, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3518, 3, 13, 14,
                                                                       1073, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3524, 3, 14, 15,
                                                                       1076, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3530, 3, 15, 16,
                                                                       1079, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3536, 3, 16, 17,
                                                                       1082, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3542, 3, 17, 18,
                                                                       1085, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3548, 3, 18, 19,
                                                                       1088, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3554, 3, 19, 20,
                                                                       1091, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3560, 0, 3, 3488,
                                                                       1058, 3494, 1094, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3578, 0, 3, 3494,
                                                                       1061, 3500, 1103, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3596, 0, 3, 3500,
                                                                       1064, 3506, 1112, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3614, 0, 3, 3506,
                                                                       1067, 3512, 1121, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3632, 0, 3, 3512,
                                                                       1070, 3518, 1130, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3650, 0, 3, 3518,
                                                                       1073, 3524, 1139, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3668, 0, 3, 3524,
                                                                       1076, 3530, 1148, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3686, 0, 3, 3530,
                                                                       1079, 3536, 1157, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3704, 0, 3, 3536,
                                                                       1082, 3542, 1166, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3722, 0, 3, 3542,
                                                                       1085, 3548, 1175, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3740, 0, 3, 3548,
                                                                       1088, 3554, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3758, 0, 3, 3560,
                                                                       1094, 3578, 61, 67, 1193,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3794, 0, 3, 3578,
                                                                       1103, 3596, 67, 73, 1211,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3830, 0, 3, 3596,
                                                                       1112, 3614, 73, 79, 1229,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3866, 0, 3, 3614,
                                                                       1121, 3632, 79, 85, 1247,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3902, 0, 3, 3632,
                                                                       1130, 3650, 85, 91, 1265,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3938, 0, 3, 3650,
                                                                       1139, 3668, 91, 97, 1283,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3974, 0, 3, 3668,
                                                                       1148, 3686, 97, 103, 1301,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4010, 0, 3, 3686,
                                                                       1157, 3704, 103, 109,
                                                                       1319, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4046, 0, 3, 3704,
                                                                       1166, 3722, 109, 115,
                                                                       1337, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4082, 0, 3, 3722,
                                                                       1175, 3740, 115, 121,
                                                                       1355, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4118, 0, 3, 3758,
                                                                       1193, 3794, 133, 143,
                                                                       1373, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4178, 0, 3, 3794,
                                                                       1211, 3830, 143, 153,
                                                                       1403, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4238, 0, 3, 3830,
                                                                       1229, 3866, 153, 163,
                                                                       1433, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4298, 0, 3, 3866,
                                                                       1247, 3902, 163, 173,
                                                                       1463, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4358, 0, 3, 3902,
                                                                       1265, 3938, 173, 183,
                                                                       1493, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4418, 0, 3, 3938,
                                                                       1283, 3974, 183, 193,
                                                                       1523, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4478, 0, 3, 3974,
                                                                       1301, 4010, 193, 203,
                                                                       1553, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4538, 0, 3, 4010,
                                                                       1319, 4046, 203, 213,
                                                                       1583, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4598, 0, 3, 4046,
                                                                       1337, 4082, 213, 223,
                                                                       1613, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4658, 0, 3, 4118,
                                                                       1373, 4178, 243, 258,
                                                                       1643, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4748, 0, 3, 4178,
                                                                       1403, 4238, 258, 273,
                                                                       1688, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4838, 0, 3, 4238,
                                                                       1433, 4298, 273, 288,
                                                                       1733, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4928, 0, 3, 4298,
                                                                       1463, 4358, 288, 303,
                                                                       1778, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5018, 0, 3, 4358,
                                                                       1493, 4418, 303, 318,
                                                                       1823, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5108, 0, 3, 4418,
                                                                       1523, 4478, 318, 333,
                                                                       1868, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5198, 0, 3, 4478,
                                                                       1553, 4538, 333, 348,
                                                                       1913, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5288, 0, 3, 4538,
                                                                       1583, 4598, 348, 363,
                                                                       1958, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5378, 0, 3, 4658,
                                                                       1643, 4748, 393, 414,
                                                                       2003, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5504, 0, 3, 4748,
                                                                       1688, 4838, 414, 435,
                                                                       2066, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5630, 0, 3, 4838,
                                                                       1733, 4928, 435, 456,
                                                                       2129, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5756, 0, 3, 4928,
                                                                       1778, 5018, 456, 477,
                                                                       2192, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5882, 0, 3, 5018,
                                                                       1823, 5108, 477, 498,
                                                                       2255, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6008, 0, 3, 5108,
                                                                       1868, 5198, 498, 519,
                                                                       2318, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6134, 0, 3, 5198,
                                                                       1913, 5288, 519, 540,
                                                                       2381, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 6260, 0, 3, 5378,
                                                                       2003, 5504, 582, 610,
                                                                       2444, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 6428, 0, 3, 5504,
                                                                       2066, 5630, 610, 638,
                                                                       2528, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 6596, 0, 3, 5630,
                                                                       2129, 5756, 638, 666,
                                                                       2612, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 6764, 0, 3, 5756,
                                                                       2192, 5882, 666, 694,
                                                                       2696, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 6932, 0, 3, 5882,
                                                                       2255, 6008, 694, 722,
                                                                       2780, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7100, 0, 3, 6008,
                                                                       2318, 6134, 722, 750,
                                                                       2864, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 7268, 0, 3, 6260,
                                                                       2444, 6428, 806, 842,
                                                                       2948, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 7484, 0, 3, 6428,
                                                                       2528, 6596, 842, 878,
                                                                       3056, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 7700, 0, 3, 6596,
                                                                       2612, 6764, 878, 914,
                                                                       3164, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 7916, 0, 3, 6764,
                                                                       2696, 6932, 914, 950,
                                                                       3272, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 8132, 0, 3, 6932,
                                                                       2780, 7100, 950, 986,
                                                                       3380, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8348, 3, 1058,
                                                                       1061, 3500, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8358, 3, 1061,
                                                                       1064, 3506, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8368, 3, 1064,
                                                                       1067, 3512, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8378, 3, 1067,
                                                                       1070, 3518, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8388, 3, 1070,
                                                                       1073, 3524, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8398, 3, 1073,
                                                                       1076, 3530, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8408, 3, 1076,
                                                                       1079, 3536, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8418, 3, 1079,
                                                                       1082, 3542, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8428, 3, 1082,
                                                                       1085, 3548, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8438, 3, 1085,
                                                                       1088, 3554, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8448, 0, 3, 8348,
                                                                       3500, 8358, 1094, 1103,
                                                                       3596, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8478, 0, 3, 8358,
                                                                       3506, 8368, 1103, 1112,
                                                                       3614, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8508, 0, 3, 8368,
                                                                       3512, 8378, 1112, 1121,
                                                                       3632, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8538, 0, 3, 8378,
                                                                       3518, 8388, 1121, 1130,
                                                                       3650, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8568, 0, 3, 8388,
                                                                       3524, 8398, 1130, 1139,
                                                                       3668, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8598, 0, 3, 8398,
                                                                       3530, 8408, 1139, 1148,
                                                                       3686, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8628, 0, 3, 8408,
                                                                       3536, 8418, 1148, 1157,
                                                                       3704, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8658, 0, 3, 8418,
                                                                       3542, 8428, 1157, 1166,
                                                                       3722, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8688, 0, 3, 8428,
                                                                       3548, 8438, 1166, 1175,
                                                                       3740, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8718, 0, 3, 8448,
                                                                       3596, 8478, 1193, 1211,
                                                                       3830, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8778, 0, 3, 8478,
                                                                       3614, 8508, 1211, 1229,
                                                                       3866, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8838, 0, 3, 8508,
                                                                       3632, 8538, 1229, 1247,
                                                                       3902, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8898, 0, 3, 8538,
                                                                       3650, 8568, 1247, 1265,
                                                                       3938, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8958, 0, 3, 8568,
                                                                       3668, 8598, 1265, 1283,
                                                                       3974, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9018, 0, 3, 8598,
                                                                       3686, 8628, 1283, 1301,
                                                                       4010, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9078, 0, 3, 8628,
                                                                       3704, 8658, 1301, 1319,
                                                                       4046, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9138, 0, 3, 8658,
                                                                       3722, 8688, 1319, 1337,
                                                                       4082, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9198, 0, 3, 8718,
                                                                       3830, 8778, 1373, 1403,
                                                                       4238, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9298, 0, 3, 8778,
                                                                       3866, 8838, 1403, 1433,
                                                                       4298, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9398, 0, 3, 8838,
                                                                       3902, 8898, 1433, 1463,
                                                                       4358, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9498, 0, 3, 8898,
                                                                       3938, 8958, 1463, 1493,
                                                                       4418, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9598, 0, 3, 8958,
                                                                       3974, 9018, 1493, 1523,
                                                                       4478, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9698, 0, 3, 9018,
                                                                       4010, 9078, 1523, 1553,
                                                                       4538, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9798, 0, 3, 9078,
                                                                       4046, 9138, 1553, 1583,
                                                                       4598, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 9898, 0, 3, 9198,
                                                                       4238, 9298, 1643, 1688,
                                                                       4838, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 10048, 0, 3, 9298,
                                                                       4298, 9398, 1688, 1733,
                                                                       4928, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 10198, 0, 3, 9398,
                                                                       4358, 9498, 1733, 1778,
                                                                       5018, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 10348, 0, 3, 9498,
                                                                       4418, 9598, 1778, 1823,
                                                                       5108, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 10498, 0, 3, 9598,
                                                                       4478, 9698, 1823, 1868,
                                                                       5198, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 10648, 0, 3, 9698,
                                                                       4538, 9798, 1868, 1913,
                                                                       5288, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 10798, 0, 3, 9898,
                                                                       4838, 10048, 2003, 2066,
                                                                       5630, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 11008, 0, 3,
                                                                       10048, 4928, 10198, 2066,
                                                                       2129, 5756, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 11218, 0, 3,
                                                                       10198, 5018, 10348, 2129,
                                                                       2192, 5882, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 11428, 0, 3,
                                                                       10348, 5108, 10498, 2192,
                                                                       2255, 6008, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 11638, 0, 3,
                                                                       10498, 5198, 10648, 2255,
                                                                       2318, 6134, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 11848, 0, 3,
                                                                       10798, 5630, 11008, 2444,
                                                                       2528, 6596, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 12128, 0, 3,
                                                                       11008, 5756, 11218, 2528,
                                                                       2612, 6764, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 12408, 0, 3,
                                                                       11218, 5882, 11428, 2612,
                                                                       2696, 6932, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 12688, 0, 3,
                                                                       11428, 6008, 11638, 2696,
                                                                       2780, 7100, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 12968, 0, 3,
                                                                       11848, 6596, 12128, 2948,
                                                                       3056, 7700, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 13328, 0, 3,
                                                                       12128, 6764, 12408, 3056,
                                                                       3164, 7916, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 13688, 0, 3,
                                                                       12408, 6932, 12688, 3164,
                                                                       3272, 8132, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14048, 3, 3488,
                                                                       3494, 8348, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14063, 3, 3494,
                                                                       3500, 8358, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14078, 3, 3500,
                                                                       3506, 8368, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14093, 3, 3506,
                                                                       3512, 8378, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14108, 3, 3512,
                                                                       3518, 8388, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14123, 3, 3518,
                                                                       3524, 8398, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14138, 3, 3524,
                                                                       3530, 8408, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14153, 3, 3530,
                                                                       3536, 8418, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14168, 3, 3536,
                                                                       3542, 8428, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14183, 3, 3542,
                                                                       3548, 8438, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14198, 0, 3,
                                                                       14048, 8348, 14063, 3560,
                                                                       3578, 8448, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14243, 0, 3,
                                                                       14063, 8358, 14078, 3578,
                                                                       3596, 8478, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14288, 0, 3,
                                                                       14078, 8368, 14093, 3596,
                                                                       3614, 8508, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14333, 0, 3,
                                                                       14093, 8378, 14108, 3614,
                                                                       3632, 8538, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14378, 0, 3,
                                                                       14108, 8388, 14123, 3632,
                                                                       3650, 8568, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14423, 0, 3,
                                                                       14123, 8398, 14138, 3650,
                                                                       3668, 8598, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14468, 0, 3,
                                                                       14138, 8408, 14153, 3668,
                                                                       3686, 8628, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14513, 0, 3,
                                                                       14153, 8418, 14168, 3686,
                                                                       3704, 8658, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14558, 0, 3,
                                                                       14168, 8428, 14183, 3704,
                                                                       3722, 8688, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 14603, 0, 3,
                                                                       14198, 8448, 14243, 3758,
                                                                       3794, 8718, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 14693, 0, 3,
                                                                       14243, 8478, 14288, 3794,
                                                                       3830, 8778, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 14783, 0, 3,
                                                                       14288, 8508, 14333, 3830,
                                                                       3866, 8838, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 14873, 0, 3,
                                                                       14333, 8538, 14378, 3866,
                                                                       3902, 8898, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 14963, 0, 3,
                                                                       14378, 8568, 14423, 3902,
                                                                       3938, 8958, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 15053, 0, 3,
                                                                       14423, 8598, 14468, 3938,
                                                                       3974, 9018, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 15143, 0, 3,
                                                                       14468, 8628, 14513, 3974,
                                                                       4010, 9078, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 15233, 0, 3,
                                                                       14513, 8658, 14558, 4010,
                                                                       4046, 9138, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 15323, 0, 3,
                                                                       14603, 8718, 14693, 4118,
                                                                       4178, 9198, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 15473, 0, 3,
                                                                       14693, 8778, 14783, 4178,
                                                                       4238, 9298, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 15623, 0, 3,
                                                                       14783, 8838, 14873, 4238,
                                                                       4298, 9398, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 15773, 0, 3,
                                                                       14873, 8898, 14963, 4298,
                                                                       4358, 9498, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 15923, 0, 3,
                                                                       14963, 8958, 15053, 4358,
                                                                       4418, 9598, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 16073, 0, 3,
                                                                       15053, 9018, 15143, 4418,
                                                                       4478, 9698, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 16223, 0, 3,
                                                                       15143, 9078, 15233, 4478,
                                                                       4538, 9798, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 16373, 0, 3,
                                                                       15323, 9198, 15473, 4658,
                                                                       4748, 9898, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 16598, 0, 3,
                                                                       15473, 9298, 15623, 4748,
                                                                       4838, 10048, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 16823, 0, 3,
                                                                       15623, 9398, 15773, 4838,
                                                                       4928, 10198, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 17048, 0, 3,
                                                                       15773, 9498, 15923, 4928,
                                                                       5018, 10348, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 17273, 0, 3,
                                                                       15923, 9598, 16073, 5018,
                                                                       5108, 10498, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 17498, 0, 3,
                                                                       16073, 9698, 16223, 5108,
                                                                       5198, 10648, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 17723, 0, 3,
                                                                       16373, 9898, 16598, 5378,
                                                                       5504, 10798, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 18038, 0, 3,
                                                                       16598, 10048, 16823, 5504,
                                                                       5630, 11008, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 18353, 0, 3,
                                                                       16823, 10198, 17048, 5630,
                                                                       5756, 11218, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 18668, 0, 3,
                                                                       17048, 10348, 17273, 5756,
                                                                       5882, 11428, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 18983, 0, 3,
                                                                       17273, 10498, 17498, 5882,
                                                                       6008, 11638, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 19298, 0, 3,
                                                                       17723, 10798, 18038, 6260,
                                                                       6428, 11848, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 19718, 0, 3,
                                                                       18038, 11008, 18353, 6428,
                                                                       6596, 12128, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 20138, 0, 3,
                                                                       18353, 11218, 18668, 6596,
                                                                       6764, 12408, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 20558, 0, 3,
                                                                       18668, 11428, 18983, 6764,
                                                                       6932, 12688, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 20978, 0, 3,
                                                                       19298, 11848, 19718, 7268,
                                                                       7484, 12968, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 21518, 0, 3,
                                                                       19718, 12128, 20138, 7484,
                                                                       7700, 13328, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 22058, 0, 3,
                                                                       20138, 12408, 20558, 7700,
                                                                       7916, 13688, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22598, 3, 8348,
                                                                       8358, 14078, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22619, 3, 8358,
                                                                       8368, 14093, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22640, 3, 8368,
                                                                       8378, 14108, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22661, 3, 8378,
                                                                       8388, 14123, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22682, 3, 8388,
                                                                       8398, 14138, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22703, 3, 8398,
                                                                       8408, 14153, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22724, 3, 8408,
                                                                       8418, 14168, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22745, 3, 8418,
                                                                       8428, 14183, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 22766, 0, 3,
                                                                       22598, 14078, 22619, 8448,
                                                                       8478, 14288, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 22829, 0, 3,
                                                                       22619, 14093, 22640, 8478,
                                                                       8508, 14333, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 22892, 0, 3,
                                                                       22640, 14108, 22661, 8508,
                                                                       8538, 14378, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 22955, 0, 3,
                                                                       22661, 14123, 22682, 8538,
                                                                       8568, 14423, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 23018, 0, 3,
                                                                       22682, 14138, 22703, 8568,
                                                                       8598, 14468, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 23081, 0, 3,
                                                                       22703, 14153, 22724, 8598,
                                                                       8628, 14513, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 23144, 0, 3,
                                                                       22724, 14168, 22745, 8628,
                                                                       8658, 14558, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 23207, 0, 3,
                                                                       22766, 14288, 22829, 8718,
                                                                       8778, 14783, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 23333, 0, 3,
                                                                       22829, 14333, 22892, 8778,
                                                                       8838, 14873, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 23459, 0, 3,
                                                                       22892, 14378, 22955, 8838,
                                                                       8898, 14963, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 23585, 0, 3,
                                                                       22955, 14423, 23018, 8898,
                                                                       8958, 15053, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 23711, 0, 3,
                                                                       23018, 14468, 23081, 8958,
                                                                       9018, 15143, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 23837, 0, 3,
                                                                       23081, 14513, 23144, 9018,
                                                                       9078, 15233, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 23963, 0, 3,
                                                                       23207, 14783, 23333, 9198,
                                                                       9298, 15623, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 24173, 0, 3,
                                                                       23333, 14873, 23459, 9298,
                                                                       9398, 15773, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 24383, 0, 3,
                                                                       23459, 14963, 23585, 9398,
                                                                       9498, 15923, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 24593, 0, 3,
                                                                       23585, 15053, 23711, 9498,
                                                                       9598, 16073, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 24803, 0, 3,
                                                                       23711, 15143, 23837, 9598,
                                                                       9698, 16223, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 25013, 0, 3,
                                                                       23963, 15623, 24173, 9898,
                                                                       10048, 16823, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 25328, 0, 3,
                                                                       24173, 15773, 24383,
                                                                       10048, 10198, 17048,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 25643, 0, 3,
                                                                       24383, 15923, 24593,
                                                                       10198, 10348, 17273,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 25958, 0, 3,
                                                                       24593, 16073, 24803,
                                                                       10348, 10498, 17498,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 26273, 0, 3,
                                                                       25013, 16823, 25328,
                                                                       10798, 11008, 18353,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 26714, 0, 3,
                                                                       25328, 17048, 25643,
                                                                       11008, 11218, 18668,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 27155, 0, 3,
                                                                       25643, 17273, 25958,
                                                                       11218, 11428, 18983,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 27596, 0, 3,
                                                                       26273, 18353, 26714,
                                                                       11848, 12128, 20138,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 28184, 0, 3,
                                                                       26714, 18668, 27155,
                                                                       12128, 12408, 20558,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 28772, 0, 3,
                                                                       27596, 20138, 28184,
                                                                       12968, 13328, 22058,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29528, 3, 14048,
                                                                       14063, 22598, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29556, 3, 14063,
                                                                       14078, 22619, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29584, 3, 14078,
                                                                       14093, 22640, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29612, 3, 14093,
                                                                       14108, 22661, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29640, 3, 14108,
                                                                       14123, 22682, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29668, 3, 14123,
                                                                       14138, 22703, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29696, 3, 14138,
                                                                       14153, 22724, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29724, 3, 14153,
                                                                       14168, 22745, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 29752, 0, 3,
                                                                       29528, 22598, 29556,
                                                                       14198, 14243, 22766,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 29836, 0, 3,
                                                                       29556, 22619, 29584,
                                                                       14243, 14288, 22829,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 29920, 0, 3,
                                                                       29584, 22640, 29612,
                                                                       14288, 14333, 22892,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 30004, 0, 3,
                                                                       29612, 22661, 29640,
                                                                       14333, 14378, 22955,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 30088, 0, 3,
                                                                       29640, 22682, 29668,
                                                                       14378, 14423, 23018,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 30172, 0, 3,
                                                                       29668, 22703, 29696,
                                                                       14423, 14468, 23081,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 30256, 0, 3,
                                                                       29696, 22724, 29724,
                                                                       14468, 14513, 23144,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 30340, 0, 3,
                                                                       29752, 22766, 29836,
                                                                       14603, 14693, 23207,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 30508, 0, 3,
                                                                       29836, 22829, 29920,
                                                                       14693, 14783, 23333,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 30676, 0, 3,
                                                                       29920, 22892, 30004,
                                                                       14783, 14873, 23459,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 30844, 0, 3,
                                                                       30004, 22955, 30088,
                                                                       14873, 14963, 23585,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 31012, 0, 3,
                                                                       30088, 23018, 30172,
                                                                       14963, 15053, 23711,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 31180, 0, 3,
                                                                       30172, 23081, 30256,
                                                                       15053, 15143, 23837,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 31348, 0, 3,
                                                                       30340, 23207, 30508,
                                                                       15323, 15473, 23963,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 31628, 0, 3,
                                                                       30508, 23333, 30676,
                                                                       15473, 15623, 24173,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 31908, 0, 3,
                                                                       30676, 23459, 30844,
                                                                       15623, 15773, 24383,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 32188, 0, 3,
                                                                       30844, 23585, 31012,
                                                                       15773, 15923, 24593,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 32468, 0, 3,
                                                                       31012, 23711, 31180,
                                                                       15923, 16073, 24803,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 32748, 0, 3,
                                                                       31348, 23963, 31628,
                                                                       16373, 16598, 25013,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 33168, 0, 3,
                                                                       31628, 24173, 31908,
                                                                       16598, 16823, 25328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 33588, 0, 3,
                                                                       31908, 24383, 32188,
                                                                       16823, 17048, 25643,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 34008, 0, 3,
                                                                       32188, 24593, 32468,
                                                                       17048, 17273, 25958,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 34428, 0, 3,
                                                                       32748, 25013, 33168,
                                                                       17723, 18038, 26273,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 35016, 0, 3,
                                                                       33168, 25328, 33588,
                                                                       18038, 18353, 26714,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 35604, 0, 3,
                                                                       33588, 25643, 34008,
                                                                       18353, 18668, 27155,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 36192, 0, 3,
                                                                       34428, 26273, 35016,
                                                                       19298, 19718, 27596,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 36976, 0, 3,
                                                                       35016, 26714, 35604,
                                                                       19718, 20138, 28184,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 37760, 0, 3,
                                                                       36192, 27596, 36976,
                                                                       20978, 21518, 28772,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 38768, 36192, 784, ncols);

                    simdfunc::contract_primitives(buffer, 39916, 37760, 1008, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 39552, 38768, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 40924, 39916, 36, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 41392, 39552, 40924, 13,
                                             nmax);

        simdtrf::transform_p_inner(buffer, 42484, 41392, 28, 13, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 42484, 39, nmax);
    }

    for (size_t m = 0; m < 507; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
