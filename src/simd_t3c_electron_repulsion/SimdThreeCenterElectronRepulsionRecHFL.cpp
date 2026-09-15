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


#include "SimdThreeCenterElectronRepulsionRecHFL.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferHD.hpp"
#include "SimdTransferHF.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_hfl_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_hfl_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 148045, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1309 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 148045, 124583, 7295, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 16,
                                                             ncols, fj, 6, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 25, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 28, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 31, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 34, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 37, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 40, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 43, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 46, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 49, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 52, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 55, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 58, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 61, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 64, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 67, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 70, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 73, 0, 3, 8, 9,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 79, 0, 3, 9, 10,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 85, 0, 3, 10, 11,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 91, 0, 3, 11, 12,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 97, 0, 3, 12, 13,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 103, 0, 3, 13, 14,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 109, 0, 3, 14, 15,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 115, 0, 3, 15, 16,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 121, 0, 3, 16, 17,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 127, 0, 3, 17, 18,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 133, 0, 3, 18, 19,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 139, 0, 3, 19, 20,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 145, 0, 3, 20, 21,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 151, 0, 3, 21, 22,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 157, 0, 3, 22, 23,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 163, 0, 3, 25, 28,
                                                                       73, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 173, 0, 3, 28, 31,
                                                                       79, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 183, 0, 3, 31, 34,
                                                                       85, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 193, 0, 3, 34, 37,
                                                                       91, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 203, 0, 3, 37, 40,
                                                                       97, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 213, 0, 3, 40, 43,
                                                                       103, 109, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 223, 0, 3, 43, 46,
                                                                       109, 115, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 233, 0, 3, 46, 49,
                                                                       115, 121, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 243, 0, 3, 49, 52,
                                                                       121, 127, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 253, 0, 3, 52, 55,
                                                                       127, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 263, 0, 3, 55, 58,
                                                                       133, 139, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 273, 0, 3, 58, 61,
                                                                       139, 145, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 283, 0, 3, 61, 64,
                                                                       145, 151, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 293, 0, 3, 64, 67,
                                                                       151, 157, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 303, 0, 3, 73, 79,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 318, 0, 3, 79, 85,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 333, 0, 3, 85, 91,
                                                                       183, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 348, 0, 3, 91, 97,
                                                                       193, 203, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 363, 0, 3, 97,
                                                                       103, 203, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 378, 0, 3, 103,
                                                                       109, 213, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 393, 0, 3, 109,
                                                                       115, 223, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 408, 0, 3, 115,
                                                                       121, 233, 243, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 423, 0, 3, 121,
                                                                       127, 243, 253, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 438, 0, 3, 127,
                                                                       133, 253, 263, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 453, 0, 3, 133,
                                                                       139, 263, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 468, 0, 3, 139,
                                                                       145, 273, 283, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 483, 0, 3, 145,
                                                                       151, 283, 293, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 498, 0, 3, 163,
                                                                       173, 303, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 519, 0, 3, 173,
                                                                       183, 318, 333, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 540, 0, 3, 183,
                                                                       193, 333, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 561, 0, 3, 193,
                                                                       203, 348, 363, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 582, 0, 3, 203,
                                                                       213, 363, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 603, 0, 3, 213,
                                                                       223, 378, 393, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 624, 0, 3, 223,
                                                                       233, 393, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 645, 0, 3, 233,
                                                                       243, 408, 423, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 666, 0, 3, 243,
                                                                       253, 423, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 687, 0, 3, 253,
                                                                       263, 438, 453, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 708, 0, 3, 263,
                                                                       273, 453, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 729, 0, 3, 273,
                                                                       283, 468, 483, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 750, 0, 3, 303,
                                                                       318, 498, 519, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 778, 0, 3, 318,
                                                                       333, 519, 540, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 806, 0, 3, 333,
                                                                       348, 540, 561, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 834, 0, 3, 348,
                                                                       363, 561, 582, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 862, 0, 3, 363,
                                                                       378, 582, 603, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 890, 0, 3, 378,
                                                                       393, 603, 624, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 918, 0, 3, 393,
                                                                       408, 624, 645, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 946, 0, 3, 408,
                                                                       423, 645, 666, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 974, 0, 3, 423,
                                                                       438, 666, 687, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1002, 0, 3, 438,
                                                                       453, 687, 708, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1030, 0, 3, 453,
                                                                       468, 708, 729, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1058, 0, 3, 498,
                                                                       519, 750, 778, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1094, 0, 3, 519,
                                                                       540, 778, 806, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1130, 0, 3, 540,
                                                                       561, 806, 834, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1166, 0, 3, 561,
                                                                       582, 834, 862, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1202, 0, 3, 582,
                                                                       603, 862, 890, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1238, 0, 3, 603,
                                                                       624, 890, 918, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1274, 0, 3, 624,
                                                                       645, 918, 946, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1310, 0, 3, 645,
                                                                       666, 946, 974, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1346, 0, 3, 666,
                                                                       687, 974, 1002, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1382, 0, 3, 687,
                                                                       708, 1002, 1030, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1418, 0, 3, 750,
                                                                       778, 1058, 1094, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1463, 0, 3, 778,
                                                                       806, 1094, 1130, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1508, 0, 3, 806,
                                                                       834, 1130, 1166, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1553, 0, 3, 834,
                                                                       862, 1166, 1202, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1598, 0, 3, 862,
                                                                       890, 1202, 1238, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1643, 0, 3, 890,
                                                                       918, 1238, 1274, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 918,
                                                                       946, 1274, 1310, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1733, 0, 3, 946,
                                                                       974, 1310, 1346, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1778, 0, 3, 974,
                                                                       1002, 1346, 1382, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1823, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1826, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1829, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1832, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1835, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1838, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1841, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1844, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1847, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1850, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1853, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1856, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1859, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1862, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1865, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1868, 3, 10, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1877, 3, 11, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1886, 3, 12, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1895, 3, 13, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1904, 3, 14, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1913, 3, 15, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1922, 3, 16, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1931, 3, 17, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1940, 3, 18, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1949, 3, 19, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1958, 3, 20, 61,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1967, 3, 21, 64,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1976, 3, 22, 67,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1985, 3, 23, 70,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1994, 3, 31, 85,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2012, 3, 34, 91,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2030, 3, 37, 97,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2048, 3, 40, 103,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2066, 3, 43, 109,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2084, 3, 46, 115,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2102, 3, 49, 121,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2120, 3, 52, 127,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2138, 3, 55, 133,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2156, 3, 58, 139,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2174, 3, 61, 145,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2192, 3, 64, 151,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2210, 3, 67, 157,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2228, 3, 85, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2258, 3, 91, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2288, 3, 97, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2318, 3, 103, 213,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2348, 3, 109, 223,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2378, 3, 115, 233,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2408, 3, 121, 243,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2438, 3, 127, 253,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2468, 3, 133, 263,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2498, 3, 139, 273,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2528, 3, 145, 283,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2558, 3, 151, 293,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2588, 3, 183, 333,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2633, 3, 193, 348,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2678, 3, 203, 363,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2723, 3, 213, 378,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2768, 3, 223, 393,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2813, 3, 233, 408,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2858, 3, 243, 423,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2903, 3, 253, 438,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2948, 3, 263, 453,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2993, 3, 273, 468,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3038, 3, 283, 483,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3083, 3, 333, 540,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3146, 3, 348, 561,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3209, 3, 363, 582,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3272, 3, 378, 603,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3335, 3, 393, 624,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3398, 3, 408, 645,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3461, 3, 423, 666,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3524, 3, 438, 687,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3587, 3, 453, 708,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3650, 3, 468, 729,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3713, 3, 540, 806,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3797, 3, 561, 834,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3881, 3, 582, 862,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3965, 3, 603, 890,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4049, 3, 624, 918,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4133, 3, 645, 946,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4217, 3, 666, 974,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4301, 3, 687,
                                                                       1002, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4385, 3, 708,
                                                                       1030, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4469, 3, 806,
                                                                       1130, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4577, 3, 834,
                                                                       1166, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4685, 3, 862,
                                                                       1202, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4793, 3, 890,
                                                                       1238, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4901, 3, 918,
                                                                       1274, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5009, 3, 946,
                                                                       1310, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5117, 3, 974,
                                                                       1346, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5225, 3, 1002,
                                                                       1382, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5333, 3, 1130,
                                                                       1508, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5468, 3, 1166,
                                                                       1553, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5603, 3, 1202,
                                                                       1598, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5738, 3, 1238,
                                                                       1643, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5873, 3, 1274,
                                                                       1688, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6008, 3, 1310,
                                                                       1733, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6143, 3, 1346,
                                                                       1778, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6278, 3, 8, 9,
                                                                       1823, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6284, 3, 9, 10,
                                                                       1826, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6290, 3, 10, 11,
                                                                       1829, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6296, 3, 11, 12,
                                                                       1832, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6302, 3, 12, 13,
                                                                       1835, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6308, 3, 13, 14,
                                                                       1838, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6314, 3, 14, 15,
                                                                       1841, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6320, 3, 15, 16,
                                                                       1844, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6326, 3, 16, 17,
                                                                       1847, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6332, 3, 17, 18,
                                                                       1850, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6338, 3, 18, 19,
                                                                       1853, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6344, 3, 19, 20,
                                                                       1856, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6350, 3, 20, 21,
                                                                       1859, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6356, 3, 21, 22,
                                                                       1862, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6362, 3, 22, 23,
                                                                       1865, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6368, 0, 3, 6278,
                                                                       1823, 6284, 1868, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6386, 0, 3, 6284,
                                                                       1826, 6290, 1877, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6404, 0, 3, 6290,
                                                                       1829, 6296, 1886, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6422, 0, 3, 6296,
                                                                       1832, 6302, 1895, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6440, 0, 3, 6302,
                                                                       1835, 6308, 1904, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6458, 0, 3, 6308,
                                                                       1838, 6314, 1913, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6476, 0, 3, 6314,
                                                                       1841, 6320, 1922, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6494, 0, 3, 6320,
                                                                       1844, 6326, 1931, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6512, 0, 3, 6326,
                                                                       1847, 6332, 1940, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6530, 0, 3, 6332,
                                                                       1850, 6338, 1949, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6548, 0, 3, 6338,
                                                                       1853, 6344, 1958, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6566, 0, 3, 6344,
                                                                       1856, 6350, 1967, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6584, 0, 3, 6350,
                                                                       1859, 6356, 1976, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6602, 0, 3, 6356,
                                                                       1862, 6362, 1985, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6620, 0, 3, 6368,
                                                                       1868, 6386, 73, 79, 1994,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6656, 0, 3, 6386,
                                                                       1877, 6404, 79, 85, 2012,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6692, 0, 3, 6404,
                                                                       1886, 6422, 85, 91, 2030,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6728, 0, 3, 6422,
                                                                       1895, 6440, 91, 97, 2048,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6764, 0, 3, 6440,
                                                                       1904, 6458, 97, 103, 2066,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6800, 0, 3, 6458,
                                                                       1913, 6476, 103, 109,
                                                                       2084, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6836, 0, 3, 6476,
                                                                       1922, 6494, 109, 115,
                                                                       2102, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6872, 0, 3, 6494,
                                                                       1931, 6512, 115, 121,
                                                                       2120, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6908, 0, 3, 6512,
                                                                       1940, 6530, 121, 127,
                                                                       2138, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6944, 0, 3, 6530,
                                                                       1949, 6548, 127, 133,
                                                                       2156, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6980, 0, 3, 6548,
                                                                       1958, 6566, 133, 139,
                                                                       2174, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7016, 0, 3, 6566,
                                                                       1967, 6584, 139, 145,
                                                                       2192, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7052, 0, 3, 6584,
                                                                       1976, 6602, 145, 151,
                                                                       2210, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7088, 0, 3, 6620,
                                                                       1994, 6656, 163, 173,
                                                                       2228, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7148, 0, 3, 6656,
                                                                       2012, 6692, 173, 183,
                                                                       2258, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7208, 0, 3, 6692,
                                                                       2030, 6728, 183, 193,
                                                                       2288, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7268, 0, 3, 6728,
                                                                       2048, 6764, 193, 203,
                                                                       2318, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7328, 0, 3, 6764,
                                                                       2066, 6800, 203, 213,
                                                                       2348, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7388, 0, 3, 6800,
                                                                       2084, 6836, 213, 223,
                                                                       2378, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7448, 0, 3, 6836,
                                                                       2102, 6872, 223, 233,
                                                                       2408, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7508, 0, 3, 6872,
                                                                       2120, 6908, 233, 243,
                                                                       2438, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7568, 0, 3, 6908,
                                                                       2138, 6944, 243, 253,
                                                                       2468, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7628, 0, 3, 6944,
                                                                       2156, 6980, 253, 263,
                                                                       2498, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7688, 0, 3, 6980,
                                                                       2174, 7016, 263, 273,
                                                                       2528, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7748, 0, 3, 7016,
                                                                       2192, 7052, 273, 283,
                                                                       2558, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7808, 0, 3, 7088,
                                                                       2228, 7148, 303, 318,
                                                                       2588, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7898, 0, 3, 7148,
                                                                       2258, 7208, 318, 333,
                                                                       2633, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7988, 0, 3, 7208,
                                                                       2288, 7268, 333, 348,
                                                                       2678, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8078, 0, 3, 7268,
                                                                       2318, 7328, 348, 363,
                                                                       2723, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8168, 0, 3, 7328,
                                                                       2348, 7388, 363, 378,
                                                                       2768, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8258, 0, 3, 7388,
                                                                       2378, 7448, 378, 393,
                                                                       2813, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8348, 0, 3, 7448,
                                                                       2408, 7508, 393, 408,
                                                                       2858, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8438, 0, 3, 7508,
                                                                       2438, 7568, 408, 423,
                                                                       2903, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8528, 0, 3, 7568,
                                                                       2468, 7628, 423, 438,
                                                                       2948, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8618, 0, 3, 7628,
                                                                       2498, 7688, 438, 453,
                                                                       2993, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 8708, 0, 3, 7688,
                                                                       2528, 7748, 453, 468,
                                                                       3038, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8798, 0, 3, 7808,
                                                                       2588, 7898, 498, 519,
                                                                       3083, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8924, 0, 3, 7898,
                                                                       2633, 7988, 519, 540,
                                                                       3146, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9050, 0, 3, 7988,
                                                                       2678, 8078, 540, 561,
                                                                       3209, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9176, 0, 3, 8078,
                                                                       2723, 8168, 561, 582,
                                                                       3272, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9302, 0, 3, 8168,
                                                                       2768, 8258, 582, 603,
                                                                       3335, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9428, 0, 3, 8258,
                                                                       2813, 8348, 603, 624,
                                                                       3398, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9554, 0, 3, 8348,
                                                                       2858, 8438, 624, 645,
                                                                       3461, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9680, 0, 3, 8438,
                                                                       2903, 8528, 645, 666,
                                                                       3524, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9806, 0, 3, 8528,
                                                                       2948, 8618, 666, 687,
                                                                       3587, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 9932, 0, 3, 8618,
                                                                       2993, 8708, 687, 708,
                                                                       3650, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10058, 0, 3, 8798,
                                                                       3083, 8924, 750, 778,
                                                                       3713, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10226, 0, 3, 8924,
                                                                       3146, 9050, 778, 806,
                                                                       3797, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10394, 0, 3, 9050,
                                                                       3209, 9176, 806, 834,
                                                                       3881, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10562, 0, 3, 9176,
                                                                       3272, 9302, 834, 862,
                                                                       3965, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10730, 0, 3, 9302,
                                                                       3335, 9428, 862, 890,
                                                                       4049, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10898, 0, 3, 9428,
                                                                       3398, 9554, 890, 918,
                                                                       4133, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11066, 0, 3, 9554,
                                                                       3461, 9680, 918, 946,
                                                                       4217, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11234, 0, 3, 9680,
                                                                       3524, 9806, 946, 974,
                                                                       4301, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11402, 0, 3, 9806,
                                                                       3587, 9932, 974, 1002,
                                                                       4385, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 11570, 0, 3,
                                                                       10058, 3713, 10226, 1058,
                                                                       1094, 4469, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 11786, 0, 3,
                                                                       10226, 3797, 10394, 1094,
                                                                       1130, 4577, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12002, 0, 3,
                                                                       10394, 3881, 10562, 1130,
                                                                       1166, 4685, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12218, 0, 3,
                                                                       10562, 3965, 10730, 1166,
                                                                       1202, 4793, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12434, 0, 3,
                                                                       10730, 4049, 10898, 1202,
                                                                       1238, 4901, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12650, 0, 3,
                                                                       10898, 4133, 11066, 1238,
                                                                       1274, 5009, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12866, 0, 3,
                                                                       11066, 4217, 11234, 1274,
                                                                       1310, 5117, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 13082, 0, 3,
                                                                       11234, 4301, 11402, 1310,
                                                                       1346, 5225, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 13298, 0, 3,
                                                                       11570, 4469, 11786, 1418,
                                                                       1463, 5333, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 13568, 0, 3,
                                                                       11786, 4577, 12002, 1463,
                                                                       1508, 5468, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 13838, 0, 3,
                                                                       12002, 4685, 12218, 1508,
                                                                       1553, 5603, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 14108, 0, 3,
                                                                       12218, 4793, 12434, 1553,
                                                                       1598, 5738, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 14378, 0, 3,
                                                                       12434, 4901, 12650, 1598,
                                                                       1643, 5873, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 14648, 0, 3,
                                                                       12650, 5009, 12866, 1643,
                                                                       1688, 6008, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 14918, 0, 3,
                                                                       12866, 5117, 13082, 1688,
                                                                       1733, 6143, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15188, 3, 1823,
                                                                       1826, 6290, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15198, 3, 1826,
                                                                       1829, 6296, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15208, 3, 1829,
                                                                       1832, 6302, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15218, 3, 1832,
                                                                       1835, 6308, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15228, 3, 1835,
                                                                       1838, 6314, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15238, 3, 1838,
                                                                       1841, 6320, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15248, 3, 1841,
                                                                       1844, 6326, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15258, 3, 1844,
                                                                       1847, 6332, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15268, 3, 1847,
                                                                       1850, 6338, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15278, 3, 1850,
                                                                       1853, 6344, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15288, 3, 1853,
                                                                       1856, 6350, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15298, 3, 1856,
                                                                       1859, 6356, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 15308, 3, 1859,
                                                                       1862, 6362, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15318, 0, 3,
                                                                       15188, 6290, 15198, 1868,
                                                                       1877, 6404, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15348, 0, 3,
                                                                       15198, 6296, 15208, 1877,
                                                                       1886, 6422, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15378, 0, 3,
                                                                       15208, 6302, 15218, 1886,
                                                                       1895, 6440, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15408, 0, 3,
                                                                       15218, 6308, 15228, 1895,
                                                                       1904, 6458, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15438, 0, 3,
                                                                       15228, 6314, 15238, 1904,
                                                                       1913, 6476, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15468, 0, 3,
                                                                       15238, 6320, 15248, 1913,
                                                                       1922, 6494, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15498, 0, 3,
                                                                       15248, 6326, 15258, 1922,
                                                                       1931, 6512, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15528, 0, 3,
                                                                       15258, 6332, 15268, 1931,
                                                                       1940, 6530, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15558, 0, 3,
                                                                       15268, 6338, 15278, 1940,
                                                                       1949, 6548, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15588, 0, 3,
                                                                       15278, 6344, 15288, 1949,
                                                                       1958, 6566, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15618, 0, 3,
                                                                       15288, 6350, 15298, 1958,
                                                                       1967, 6584, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15648, 0, 3,
                                                                       15298, 6356, 15308, 1967,
                                                                       1976, 6602, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 15678, 0, 3,
                                                                       15318, 6404, 15348, 1994,
                                                                       2012, 6692, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 15738, 0, 3,
                                                                       15348, 6422, 15378, 2012,
                                                                       2030, 6728, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 15798, 0, 3,
                                                                       15378, 6440, 15408, 2030,
                                                                       2048, 6764, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 15858, 0, 3,
                                                                       15408, 6458, 15438, 2048,
                                                                       2066, 6800, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 15918, 0, 3,
                                                                       15438, 6476, 15468, 2066,
                                                                       2084, 6836, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 15978, 0, 3,
                                                                       15468, 6494, 15498, 2084,
                                                                       2102, 6872, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16038, 0, 3,
                                                                       15498, 6512, 15528, 2102,
                                                                       2120, 6908, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16098, 0, 3,
                                                                       15528, 6530, 15558, 2120,
                                                                       2138, 6944, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16158, 0, 3,
                                                                       15558, 6548, 15588, 2138,
                                                                       2156, 6980, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16218, 0, 3,
                                                                       15588, 6566, 15618, 2156,
                                                                       2174, 7016, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16278, 0, 3,
                                                                       15618, 6584, 15648, 2174,
                                                                       2192, 7052, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 16338, 0, 3,
                                                                       15678, 6692, 15738, 2228,
                                                                       2258, 7208, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 16438, 0, 3,
                                                                       15738, 6728, 15798, 2258,
                                                                       2288, 7268, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 16538, 0, 3,
                                                                       15798, 6764, 15858, 2288,
                                                                       2318, 7328, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 16638, 0, 3,
                                                                       15858, 6800, 15918, 2318,
                                                                       2348, 7388, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 16738, 0, 3,
                                                                       15918, 6836, 15978, 2348,
                                                                       2378, 7448, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 16838, 0, 3,
                                                                       15978, 6872, 16038, 2378,
                                                                       2408, 7508, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 16938, 0, 3,
                                                                       16038, 6908, 16098, 2408,
                                                                       2438, 7568, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17038, 0, 3,
                                                                       16098, 6944, 16158, 2438,
                                                                       2468, 7628, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17138, 0, 3,
                                                                       16158, 6980, 16218, 2468,
                                                                       2498, 7688, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 17238, 0, 3,
                                                                       16218, 7016, 16278, 2498,
                                                                       2528, 7748, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17338, 0, 3,
                                                                       16338, 7208, 16438, 2588,
                                                                       2633, 7988, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17488, 0, 3,
                                                                       16438, 7268, 16538, 2633,
                                                                       2678, 8078, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17638, 0, 3,
                                                                       16538, 7328, 16638, 2678,
                                                                       2723, 8168, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17788, 0, 3,
                                                                       16638, 7388, 16738, 2723,
                                                                       2768, 8258, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17938, 0, 3,
                                                                       16738, 7448, 16838, 2768,
                                                                       2813, 8348, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18088, 0, 3,
                                                                       16838, 7508, 16938, 2813,
                                                                       2858, 8438, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18238, 0, 3,
                                                                       16938, 7568, 17038, 2858,
                                                                       2903, 8528, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18388, 0, 3,
                                                                       17038, 7628, 17138, 2903,
                                                                       2948, 8618, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 18538, 0, 3,
                                                                       17138, 7688, 17238, 2948,
                                                                       2993, 8708, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 18688, 0, 3,
                                                                       17338, 7988, 17488, 3083,
                                                                       3146, 9050, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 18898, 0, 3,
                                                                       17488, 8078, 17638, 3146,
                                                                       3209, 9176, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19108, 0, 3,
                                                                       17638, 8168, 17788, 3209,
                                                                       3272, 9302, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19318, 0, 3,
                                                                       17788, 8258, 17938, 3272,
                                                                       3335, 9428, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19528, 0, 3,
                                                                       17938, 8348, 18088, 3335,
                                                                       3398, 9554, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19738, 0, 3,
                                                                       18088, 8438, 18238, 3398,
                                                                       3461, 9680, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 19948, 0, 3,
                                                                       18238, 8528, 18388, 3461,
                                                                       3524, 9806, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 20158, 0, 3,
                                                                       18388, 8618, 18538, 3524,
                                                                       3587, 9932, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 20368, 0, 3,
                                                                       18688, 9050, 18898, 3713,
                                                                       3797, 10394, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 20648, 0, 3,
                                                                       18898, 9176, 19108, 3797,
                                                                       3881, 10562, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 20928, 0, 3,
                                                                       19108, 9302, 19318, 3881,
                                                                       3965, 10730, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 21208, 0, 3,
                                                                       19318, 9428, 19528, 3965,
                                                                       4049, 10898, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 21488, 0, 3,
                                                                       19528, 9554, 19738, 4049,
                                                                       4133, 11066, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 21768, 0, 3,
                                                                       19738, 9680, 19948, 4133,
                                                                       4217, 11234, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 22048, 0, 3,
                                                                       19948, 9806, 20158, 4217,
                                                                       4301, 11402, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 22328, 0, 3,
                                                                       20368, 10394, 20648, 4469,
                                                                       4577, 12002, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 22688, 0, 3,
                                                                       20648, 10562, 20928, 4577,
                                                                       4685, 12218, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 23048, 0, 3,
                                                                       20928, 10730, 21208, 4685,
                                                                       4793, 12434, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 23408, 0, 3,
                                                                       21208, 10898, 21488, 4793,
                                                                       4901, 12650, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 23768, 0, 3,
                                                                       21488, 11066, 21768, 4901,
                                                                       5009, 12866, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 24128, 0, 3,
                                                                       21768, 11234, 22048, 5009,
                                                                       5117, 13082, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 24488, 0, 3,
                                                                       22328, 12002, 22688, 5333,
                                                                       5468, 13838, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 24938, 0, 3,
                                                                       22688, 12218, 23048, 5468,
                                                                       5603, 14108, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 25388, 0, 3,
                                                                       23048, 12434, 23408, 5603,
                                                                       5738, 14378, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 25838, 0, 3,
                                                                       23408, 12650, 23768, 5738,
                                                                       5873, 14648, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 26288, 0, 3,
                                                                       23768, 12866, 24128, 5873,
                                                                       6008, 14918, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26738, 3, 6278,
                                                                       6284, 15188, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26753, 3, 6284,
                                                                       6290, 15198, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26768, 3, 6290,
                                                                       6296, 15208, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26783, 3, 6296,
                                                                       6302, 15218, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26798, 3, 6302,
                                                                       6308, 15228, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26813, 3, 6308,
                                                                       6314, 15238, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26828, 3, 6314,
                                                                       6320, 15248, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26843, 3, 6320,
                                                                       6326, 15258, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26858, 3, 6326,
                                                                       6332, 15268, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26873, 3, 6332,
                                                                       6338, 15278, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26888, 3, 6338,
                                                                       6344, 15288, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26903, 3, 6344,
                                                                       6350, 15298, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26918, 3, 6350,
                                                                       6356, 15308, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26933, 0, 3,
                                                                       26738, 15188, 26753, 6368,
                                                                       6386, 15318, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26978, 0, 3,
                                                                       26753, 15198, 26768, 6386,
                                                                       6404, 15348, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27023, 0, 3,
                                                                       26768, 15208, 26783, 6404,
                                                                       6422, 15378, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27068, 0, 3,
                                                                       26783, 15218, 26798, 6422,
                                                                       6440, 15408, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27113, 0, 3,
                                                                       26798, 15228, 26813, 6440,
                                                                       6458, 15438, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27158, 0, 3,
                                                                       26813, 15238, 26828, 6458,
                                                                       6476, 15468, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27203, 0, 3,
                                                                       26828, 15248, 26843, 6476,
                                                                       6494, 15498, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27248, 0, 3,
                                                                       26843, 15258, 26858, 6494,
                                                                       6512, 15528, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27293, 0, 3,
                                                                       26858, 15268, 26873, 6512,
                                                                       6530, 15558, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27338, 0, 3,
                                                                       26873, 15278, 26888, 6530,
                                                                       6548, 15588, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27383, 0, 3,
                                                                       26888, 15288, 26903, 6548,
                                                                       6566, 15618, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27428, 0, 3,
                                                                       26903, 15298, 26918, 6566,
                                                                       6584, 15648, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 27473, 0, 3,
                                                                       26933, 15318, 26978, 6620,
                                                                       6656, 15678, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 27563, 0, 3,
                                                                       26978, 15348, 27023, 6656,
                                                                       6692, 15738, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 27653, 0, 3,
                                                                       27023, 15378, 27068, 6692,
                                                                       6728, 15798, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 27743, 0, 3,
                                                                       27068, 15408, 27113, 6728,
                                                                       6764, 15858, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 27833, 0, 3,
                                                                       27113, 15438, 27158, 6764,
                                                                       6800, 15918, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 27923, 0, 3,
                                                                       27158, 15468, 27203, 6800,
                                                                       6836, 15978, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 28013, 0, 3,
                                                                       27203, 15498, 27248, 6836,
                                                                       6872, 16038, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 28103, 0, 3,
                                                                       27248, 15528, 27293, 6872,
                                                                       6908, 16098, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 28193, 0, 3,
                                                                       27293, 15558, 27338, 6908,
                                                                       6944, 16158, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 28283, 0, 3,
                                                                       27338, 15588, 27383, 6944,
                                                                       6980, 16218, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 28373, 0, 3,
                                                                       27383, 15618, 27428, 6980,
                                                                       7016, 16278, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 28463, 0, 3,
                                                                       27473, 15678, 27563, 7088,
                                                                       7148, 16338, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 28613, 0, 3,
                                                                       27563, 15738, 27653, 7148,
                                                                       7208, 16438, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 28763, 0, 3,
                                                                       27653, 15798, 27743, 7208,
                                                                       7268, 16538, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 28913, 0, 3,
                                                                       27743, 15858, 27833, 7268,
                                                                       7328, 16638, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 29063, 0, 3,
                                                                       27833, 15918, 27923, 7328,
                                                                       7388, 16738, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 29213, 0, 3,
                                                                       27923, 15978, 28013, 7388,
                                                                       7448, 16838, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 29363, 0, 3,
                                                                       28013, 16038, 28103, 7448,
                                                                       7508, 16938, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 29513, 0, 3,
                                                                       28103, 16098, 28193, 7508,
                                                                       7568, 17038, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 29663, 0, 3,
                                                                       28193, 16158, 28283, 7568,
                                                                       7628, 17138, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 29813, 0, 3,
                                                                       28283, 16218, 28373, 7628,
                                                                       7688, 17238, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 29963, 0, 3,
                                                                       28463, 16338, 28613, 7808,
                                                                       7898, 17338, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 30188, 0, 3,
                                                                       28613, 16438, 28763, 7898,
                                                                       7988, 17488, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 30413, 0, 3,
                                                                       28763, 16538, 28913, 7988,
                                                                       8078, 17638, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 30638, 0, 3,
                                                                       28913, 16638, 29063, 8078,
                                                                       8168, 17788, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 30863, 0, 3,
                                                                       29063, 16738, 29213, 8168,
                                                                       8258, 17938, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 31088, 0, 3,
                                                                       29213, 16838, 29363, 8258,
                                                                       8348, 18088, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 31313, 0, 3,
                                                                       29363, 16938, 29513, 8348,
                                                                       8438, 18238, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 31538, 0, 3,
                                                                       29513, 17038, 29663, 8438,
                                                                       8528, 18388, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 31763, 0, 3,
                                                                       29663, 17138, 29813, 8528,
                                                                       8618, 18538, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 31988, 0, 3,
                                                                       29963, 17338, 30188, 8798,
                                                                       8924, 18688, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 32303, 0, 3,
                                                                       30188, 17488, 30413, 8924,
                                                                       9050, 18898, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 32618, 0, 3,
                                                                       30413, 17638, 30638, 9050,
                                                                       9176, 19108, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 32933, 0, 3,
                                                                       30638, 17788, 30863, 9176,
                                                                       9302, 19318, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 33248, 0, 3,
                                                                       30863, 17938, 31088, 9302,
                                                                       9428, 19528, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 33563, 0, 3,
                                                                       31088, 18088, 31313, 9428,
                                                                       9554, 19738, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 33878, 0, 3,
                                                                       31313, 18238, 31538, 9554,
                                                                       9680, 19948, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 34193, 0, 3,
                                                                       31538, 18388, 31763, 9680,
                                                                       9806, 20158, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 34508, 0, 3,
                                                                       31988, 18688, 32303,
                                                                       10058, 10226, 20368,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 34928, 0, 3,
                                                                       32303, 18898, 32618,
                                                                       10226, 10394, 20648,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 35348, 0, 3,
                                                                       32618, 19108, 32933,
                                                                       10394, 10562, 20928,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 35768, 0, 3,
                                                                       32933, 19318, 33248,
                                                                       10562, 10730, 21208,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 36188, 0, 3,
                                                                       33248, 19528, 33563,
                                                                       10730, 10898, 21488,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 36608, 0, 3,
                                                                       33563, 19738, 33878,
                                                                       10898, 11066, 21768,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 37028, 0, 3,
                                                                       33878, 19948, 34193,
                                                                       11066, 11234, 22048,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 37448, 0, 3,
                                                                       34508, 20368, 34928,
                                                                       11570, 11786, 22328,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 37988, 0, 3,
                                                                       34928, 20648, 35348,
                                                                       11786, 12002, 22688,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 38528, 0, 3,
                                                                       35348, 20928, 35768,
                                                                       12002, 12218, 23048,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 39068, 0, 3,
                                                                       35768, 21208, 36188,
                                                                       12218, 12434, 23408,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 39608, 0, 3,
                                                                       36188, 21488, 36608,
                                                                       12434, 12650, 23768,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 40148, 0, 3,
                                                                       36608, 21768, 37028,
                                                                       12650, 12866, 24128,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 40688, 0, 3,
                                                                       37448, 22328, 37988,
                                                                       13298, 13568, 24488,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 41363, 0, 3,
                                                                       37988, 22688, 38528,
                                                                       13568, 13838, 24938,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 42038, 0, 3,
                                                                       38528, 23048, 39068,
                                                                       13838, 14108, 25388,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 42713, 0, 3,
                                                                       39068, 23408, 39608,
                                                                       14108, 14378, 25838,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 43388, 0, 3,
                                                                       39608, 23768, 40148,
                                                                       14378, 14648, 26288,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44063, 3, 15188,
                                                                       15198, 26768, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44084, 3, 15198,
                                                                       15208, 26783, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44105, 3, 15208,
                                                                       15218, 26798, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44126, 3, 15218,
                                                                       15228, 26813, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44147, 3, 15228,
                                                                       15238, 26828, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44168, 3, 15238,
                                                                       15248, 26843, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44189, 3, 15248,
                                                                       15258, 26858, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44210, 3, 15258,
                                                                       15268, 26873, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44231, 3, 15268,
                                                                       15278, 26888, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44252, 3, 15278,
                                                                       15288, 26903, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 44273, 3, 15288,
                                                                       15298, 26918, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 44294, 0, 3,
                                                                       44063, 26768, 44084,
                                                                       15318, 15348, 27023,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 44357, 0, 3,
                                                                       44084, 26783, 44105,
                                                                       15348, 15378, 27068,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 44420, 0, 3,
                                                                       44105, 26798, 44126,
                                                                       15378, 15408, 27113,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 44483, 0, 3,
                                                                       44126, 26813, 44147,
                                                                       15408, 15438, 27158,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 44546, 0, 3,
                                                                       44147, 26828, 44168,
                                                                       15438, 15468, 27203,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 44609, 0, 3,
                                                                       44168, 26843, 44189,
                                                                       15468, 15498, 27248,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 44672, 0, 3,
                                                                       44189, 26858, 44210,
                                                                       15498, 15528, 27293,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 44735, 0, 3,
                                                                       44210, 26873, 44231,
                                                                       15528, 15558, 27338,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 44798, 0, 3,
                                                                       44231, 26888, 44252,
                                                                       15558, 15588, 27383,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 44861, 0, 3,
                                                                       44252, 26903, 44273,
                                                                       15588, 15618, 27428,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 44924, 0, 3,
                                                                       44294, 27023, 44357,
                                                                       15678, 15738, 27653,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 45050, 0, 3,
                                                                       44357, 27068, 44420,
                                                                       15738, 15798, 27743,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 45176, 0, 3,
                                                                       44420, 27113, 44483,
                                                                       15798, 15858, 27833,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 45302, 0, 3,
                                                                       44483, 27158, 44546,
                                                                       15858, 15918, 27923,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 45428, 0, 3,
                                                                       44546, 27203, 44609,
                                                                       15918, 15978, 28013,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 45554, 0, 3,
                                                                       44609, 27248, 44672,
                                                                       15978, 16038, 28103,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 45680, 0, 3,
                                                                       44672, 27293, 44735,
                                                                       16038, 16098, 28193,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 45806, 0, 3,
                                                                       44735, 27338, 44798,
                                                                       16098, 16158, 28283,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 45932, 0, 3,
                                                                       44798, 27383, 44861,
                                                                       16158, 16218, 28373,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 46058, 0, 3,
                                                                       44924, 27653, 45050,
                                                                       16338, 16438, 28763,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 46268, 0, 3,
                                                                       45050, 27743, 45176,
                                                                       16438, 16538, 28913,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 46478, 0, 3,
                                                                       45176, 27833, 45302,
                                                                       16538, 16638, 29063,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 46688, 0, 3,
                                                                       45302, 27923, 45428,
                                                                       16638, 16738, 29213,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 46898, 0, 3,
                                                                       45428, 28013, 45554,
                                                                       16738, 16838, 29363,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 47108, 0, 3,
                                                                       45554, 28103, 45680,
                                                                       16838, 16938, 29513,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 47318, 0, 3,
                                                                       45680, 28193, 45806,
                                                                       16938, 17038, 29663,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 47528, 0, 3,
                                                                       45806, 28283, 45932,
                                                                       17038, 17138, 29813,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 47738, 0, 3,
                                                                       46058, 28763, 46268,
                                                                       17338, 17488, 30413,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 48053, 0, 3,
                                                                       46268, 28913, 46478,
                                                                       17488, 17638, 30638,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 48368, 0, 3,
                                                                       46478, 29063, 46688,
                                                                       17638, 17788, 30863,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 48683, 0, 3,
                                                                       46688, 29213, 46898,
                                                                       17788, 17938, 31088,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 48998, 0, 3,
                                                                       46898, 29363, 47108,
                                                                       17938, 18088, 31313,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 49313, 0, 3,
                                                                       47108, 29513, 47318,
                                                                       18088, 18238, 31538,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 49628, 0, 3,
                                                                       47318, 29663, 47528,
                                                                       18238, 18388, 31763,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 49943, 0, 3,
                                                                       47738, 30413, 48053,
                                                                       18688, 18898, 32618,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 50384, 0, 3,
                                                                       48053, 30638, 48368,
                                                                       18898, 19108, 32933,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 50825, 0, 3,
                                                                       48368, 30863, 48683,
                                                                       19108, 19318, 33248,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 51266, 0, 3,
                                                                       48683, 31088, 48998,
                                                                       19318, 19528, 33563,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 51707, 0, 3,
                                                                       48998, 31313, 49313,
                                                                       19528, 19738, 33878,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 52148, 0, 3,
                                                                       49313, 31538, 49628,
                                                                       19738, 19948, 34193,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 52589, 0, 3,
                                                                       49943, 32618, 50384,
                                                                       20368, 20648, 35348,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 53177, 0, 3,
                                                                       50384, 32933, 50825,
                                                                       20648, 20928, 35768,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 53765, 0, 3,
                                                                       50825, 33248, 51266,
                                                                       20928, 21208, 36188,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 54353, 0, 3,
                                                                       51266, 33563, 51707,
                                                                       21208, 21488, 36608,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 54941, 0, 3,
                                                                       51707, 33878, 52148,
                                                                       21488, 21768, 37028,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 55529, 0, 3,
                                                                       52589, 35348, 53177,
                                                                       22328, 22688, 38528,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 56285, 0, 3,
                                                                       53177, 35768, 53765,
                                                                       22688, 23048, 39068,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 57041, 0, 3,
                                                                       53765, 36188, 54353,
                                                                       23048, 23408, 39608,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 57797, 0, 3,
                                                                       54353, 36608, 54941,
                                                                       23408, 23768, 40148,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 58553, 0, 3,
                                                                       55529, 38528, 56285,
                                                                       24488, 24938, 42038,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 59498, 0, 3,
                                                                       56285, 39068, 57041,
                                                                       24938, 25388, 42713,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 60443, 0, 3,
                                                                       57041, 39608, 57797,
                                                                       25388, 25838, 43388,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61388, 3, 26738,
                                                                       26753, 44063, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61416, 3, 26753,
                                                                       26768, 44084, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61444, 3, 26768,
                                                                       26783, 44105, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61472, 3, 26783,
                                                                       26798, 44126, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61500, 3, 26798,
                                                                       26813, 44147, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61528, 3, 26813,
                                                                       26828, 44168, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61556, 3, 26828,
                                                                       26843, 44189, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61584, 3, 26843,
                                                                       26858, 44210, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61612, 3, 26858,
                                                                       26873, 44231, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61640, 3, 26873,
                                                                       26888, 44252, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61668, 3, 26888,
                                                                       26903, 44273, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 61696, 0, 3,
                                                                       61388, 44063, 61416,
                                                                       26933, 26978, 44294,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 61780, 0, 3,
                                                                       61416, 44084, 61444,
                                                                       26978, 27023, 44357,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 61864, 0, 3,
                                                                       61444, 44105, 61472,
                                                                       27023, 27068, 44420,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 61948, 0, 3,
                                                                       61472, 44126, 61500,
                                                                       27068, 27113, 44483,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62032, 0, 3,
                                                                       61500, 44147, 61528,
                                                                       27113, 27158, 44546,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62116, 0, 3,
                                                                       61528, 44168, 61556,
                                                                       27158, 27203, 44609,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62200, 0, 3,
                                                                       61556, 44189, 61584,
                                                                       27203, 27248, 44672,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62284, 0, 3,
                                                                       61584, 44210, 61612,
                                                                       27248, 27293, 44735,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62368, 0, 3,
                                                                       61612, 44231, 61640,
                                                                       27293, 27338, 44798,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62452, 0, 3,
                                                                       61640, 44252, 61668,
                                                                       27338, 27383, 44861,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 62536, 0, 3,
                                                                       61696, 44294, 61780,
                                                                       27473, 27563, 44924,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 62704, 0, 3,
                                                                       61780, 44357, 61864,
                                                                       27563, 27653, 45050,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 62872, 0, 3,
                                                                       61864, 44420, 61948,
                                                                       27653, 27743, 45176,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 63040, 0, 3,
                                                                       61948, 44483, 62032,
                                                                       27743, 27833, 45302,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 63208, 0, 3,
                                                                       62032, 44546, 62116,
                                                                       27833, 27923, 45428,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 63376, 0, 3,
                                                                       62116, 44609, 62200,
                                                                       27923, 28013, 45554,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 63544, 0, 3,
                                                                       62200, 44672, 62284,
                                                                       28013, 28103, 45680,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 63712, 0, 3,
                                                                       62284, 44735, 62368,
                                                                       28103, 28193, 45806,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 63880, 0, 3,
                                                                       62368, 44798, 62452,
                                                                       28193, 28283, 45932,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 64048, 0, 3,
                                                                       62536, 44924, 62704,
                                                                       28463, 28613, 46058,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 64328, 0, 3,
                                                                       62704, 45050, 62872,
                                                                       28613, 28763, 46268,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 64608, 0, 3,
                                                                       62872, 45176, 63040,
                                                                       28763, 28913, 46478,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 64888, 0, 3,
                                                                       63040, 45302, 63208,
                                                                       28913, 29063, 46688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 65168, 0, 3,
                                                                       63208, 45428, 63376,
                                                                       29063, 29213, 46898,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 65448, 0, 3,
                                                                       63376, 45554, 63544,
                                                                       29213, 29363, 47108,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 65728, 0, 3,
                                                                       63544, 45680, 63712,
                                                                       29363, 29513, 47318,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 66008, 0, 3,
                                                                       63712, 45806, 63880,
                                                                       29513, 29663, 47528,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 66288, 0, 3,
                                                                       64048, 46058, 64328,
                                                                       29963, 30188, 47738,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 66708, 0, 3,
                                                                       64328, 46268, 64608,
                                                                       30188, 30413, 48053,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 67128, 0, 3,
                                                                       64608, 46478, 64888,
                                                                       30413, 30638, 48368,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 67548, 0, 3,
                                                                       64888, 46688, 65168,
                                                                       30638, 30863, 48683,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 67968, 0, 3,
                                                                       65168, 46898, 65448,
                                                                       30863, 31088, 48998,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 68388, 0, 3,
                                                                       65448, 47108, 65728,
                                                                       31088, 31313, 49313,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 68808, 0, 3,
                                                                       65728, 47318, 66008,
                                                                       31313, 31538, 49628,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 69228, 0, 3,
                                                                       66288, 47738, 66708,
                                                                       31988, 32303, 49943,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 69816, 0, 3,
                                                                       66708, 48053, 67128,
                                                                       32303, 32618, 50384,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 70404, 0, 3,
                                                                       67128, 48368, 67548,
                                                                       32618, 32933, 50825,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 70992, 0, 3,
                                                                       67548, 48683, 67968,
                                                                       32933, 33248, 51266,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 71580, 0, 3,
                                                                       67968, 48998, 68388,
                                                                       33248, 33563, 51707,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 72168, 0, 3,
                                                                       68388, 49313, 68808,
                                                                       33563, 33878, 52148,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 72756, 0, 3,
                                                                       69228, 49943, 69816,
                                                                       34508, 34928, 52589,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 73540, 0, 3,
                                                                       69816, 50384, 70404,
                                                                       34928, 35348, 53177,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 74324, 0, 3,
                                                                       70404, 50825, 70992,
                                                                       35348, 35768, 53765,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 75108, 0, 3,
                                                                       70992, 51266, 71580,
                                                                       35768, 36188, 54353,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 75892, 0, 3,
                                                                       71580, 51707, 72168,
                                                                       36188, 36608, 54941,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 76676, 0, 3,
                                                                       72756, 52589, 73540,
                                                                       37448, 37988, 55529,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 77684, 0, 3,
                                                                       73540, 53177, 74324,
                                                                       37988, 38528, 56285,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 78692, 0, 3,
                                                                       74324, 53765, 75108,
                                                                       38528, 39068, 57041,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 79700, 0, 3,
                                                                       75108, 54353, 75892,
                                                                       39068, 39608, 57797,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 80708, 0, 3,
                                                                       76676, 55529, 77684,
                                                                       40688, 41363, 58553,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 81968, 0, 3,
                                                                       77684, 56285, 78692,
                                                                       41363, 42038, 59498,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 83228, 0, 3,
                                                                       78692, 57041, 79700,
                                                                       42038, 42713, 60443,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 84488, 3, 44063,
                                                                       44084, 61444, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 84524, 3, 44084,
                                                                       44105, 61472, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 84560, 3, 44105,
                                                                       44126, 61500, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 84596, 3, 44126,
                                                                       44147, 61528, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 84632, 3, 44147,
                                                                       44168, 61556, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 84668, 3, 44168,
                                                                       44189, 61584, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 84704, 3, 44189,
                                                                       44210, 61612, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 84740, 3, 44210,
                                                                       44231, 61640, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 84776, 3, 44231,
                                                                       44252, 61668, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 84812, 0, 3,
                                                                       84488, 61444, 84524,
                                                                       44294, 44357, 61864,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 84920, 0, 3,
                                                                       84524, 61472, 84560,
                                                                       44357, 44420, 61948,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 85028, 0, 3,
                                                                       84560, 61500, 84596,
                                                                       44420, 44483, 62032,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 85136, 0, 3,
                                                                       84596, 61528, 84632,
                                                                       44483, 44546, 62116,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 85244, 0, 3,
                                                                       84632, 61556, 84668,
                                                                       44546, 44609, 62200,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 85352, 0, 3,
                                                                       84668, 61584, 84704,
                                                                       44609, 44672, 62284,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 85460, 0, 3,
                                                                       84704, 61612, 84740,
                                                                       44672, 44735, 62368,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 85568, 0, 3,
                                                                       84740, 61640, 84776,
                                                                       44735, 44798, 62452,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 85676, 0, 3,
                                                                       84812, 61864, 84920,
                                                                       44924, 45050, 62872,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 85892, 0, 3,
                                                                       84920, 61948, 85028,
                                                                       45050, 45176, 63040,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 86108, 0, 3,
                                                                       85028, 62032, 85136,
                                                                       45176, 45302, 63208,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 86324, 0, 3,
                                                                       85136, 62116, 85244,
                                                                       45302, 45428, 63376,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 86540, 0, 3,
                                                                       85244, 62200, 85352,
                                                                       45428, 45554, 63544,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 86756, 0, 3,
                                                                       85352, 62284, 85460,
                                                                       45554, 45680, 63712,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 86972, 0, 3,
                                                                       85460, 62368, 85568,
                                                                       45680, 45806, 63880,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 87188, 0, 3,
                                                                       85676, 62872, 85892,
                                                                       46058, 46268, 64608,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 87548, 0, 3,
                                                                       85892, 63040, 86108,
                                                                       46268, 46478, 64888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 87908, 0, 3,
                                                                       86108, 63208, 86324,
                                                                       46478, 46688, 65168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 88268, 0, 3,
                                                                       86324, 63376, 86540,
                                                                       46688, 46898, 65448,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 88628, 0, 3,
                                                                       86540, 63544, 86756,
                                                                       46898, 47108, 65728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 88988, 0, 3,
                                                                       86756, 63712, 86972,
                                                                       47108, 47318, 66008,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 89348, 0, 3,
                                                                       87188, 64608, 87548,
                                                                       47738, 48053, 67128,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 89888, 0, 3,
                                                                       87548, 64888, 87908,
                                                                       48053, 48368, 67548,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 90428, 0, 3,
                                                                       87908, 65168, 88268,
                                                                       48368, 48683, 67968,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 90968, 0, 3,
                                                                       88268, 65448, 88628,
                                                                       48683, 48998, 68388,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 91508, 0, 3,
                                                                       88628, 65728, 88988,
                                                                       48998, 49313, 68808,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 92048, 0, 3,
                                                                       89348, 67128, 89888,
                                                                       49943, 50384, 70404,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 92804, 0, 3,
                                                                       89888, 67548, 90428,
                                                                       50384, 50825, 70992,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 93560, 0, 3,
                                                                       90428, 67968, 90968,
                                                                       50825, 51266, 71580,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 94316, 0, 3,
                                                                       90968, 68388, 91508,
                                                                       51266, 51707, 72168,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 95072, 0, 3,
                                                                       92048, 70404, 92804,
                                                                       52589, 53177, 74324,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 96080, 0, 3,
                                                                       92804, 70992, 93560,
                                                                       53177, 53765, 75108,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 97088, 0, 3,
                                                                       93560, 71580, 94316,
                                                                       53765, 54353, 75892,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 98096, 0, 3,
                                                                       95072, 74324, 96080,
                                                                       55529, 56285, 78692,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 99392, 0, 3,
                                                                       96080, 75108, 97088,
                                                                       56285, 57041, 79700,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 100688, 0, 3,
                                                                       98096, 78692, 99392,
                                                                       58553, 59498, 83228,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 102308, 3, 61388,
                                                                       61416, 84488, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 102353, 3, 61416,
                                                                       61444, 84524, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 102398, 3, 61444,
                                                                       61472, 84560, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 102443, 3, 61472,
                                                                       61500, 84596, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 102488, 3, 61500,
                                                                       61528, 84632, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 102533, 3, 61528,
                                                                       61556, 84668, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 102578, 3, 61556,
                                                                       61584, 84704, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 102623, 3, 61584,
                                                                       61612, 84740, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 102668, 3, 61612,
                                                                       61640, 84776, ncols,
                                                                       gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 102713, 0, 3,
                                                                       102308, 84488, 102353,
                                                                       61696, 61780, 84812,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 102848, 0, 3,
                                                                       102353, 84524, 102398,
                                                                       61780, 61864, 84920,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 102983, 0, 3,
                                                                       102398, 84560, 102443,
                                                                       61864, 61948, 85028,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 103118, 0, 3,
                                                                       102443, 84596, 102488,
                                                                       61948, 62032, 85136,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 103253, 0, 3,
                                                                       102488, 84632, 102533,
                                                                       62032, 62116, 85244,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 103388, 0, 3,
                                                                       102533, 84668, 102578,
                                                                       62116, 62200, 85352,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 103523, 0, 3,
                                                                       102578, 84704, 102623,
                                                                       62200, 62284, 85460,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 103658, 0, 3,
                                                                       102623, 84740, 102668,
                                                                       62284, 62368, 85568,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 103793, 0, 3,
                                                                       102713, 84812, 102848,
                                                                       62536, 62704, 85676,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 104063, 0, 3,
                                                                       102848, 84920, 102983,
                                                                       62704, 62872, 85892,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 104333, 0, 3,
                                                                       102983, 85028, 103118,
                                                                       62872, 63040, 86108,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 104603, 0, 3,
                                                                       103118, 85136, 103253,
                                                                       63040, 63208, 86324,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 104873, 0, 3,
                                                                       103253, 85244, 103388,
                                                                       63208, 63376, 86540,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 105143, 0, 3,
                                                                       103388, 85352, 103523,
                                                                       63376, 63544, 86756,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 105413, 0, 3,
                                                                       103523, 85460, 103658,
                                                                       63544, 63712, 86972,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 105683, 0, 3,
                                                                       103793, 85676, 104063,
                                                                       64048, 64328, 87188,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 106133, 0, 3,
                                                                       104063, 85892, 104333,
                                                                       64328, 64608, 87548,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 106583, 0, 3,
                                                                       104333, 86108, 104603,
                                                                       64608, 64888, 87908,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 107033, 0, 3,
                                                                       104603, 86324, 104873,
                                                                       64888, 65168, 88268,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 107483, 0, 3,
                                                                       104873, 86540, 105143,
                                                                       65168, 65448, 88628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 107933, 0, 3,
                                                                       105143, 86756, 105413,
                                                                       65448, 65728, 88988,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 108383, 0, 3,
                                                                       105683, 87188, 106133,
                                                                       66288, 66708, 89348,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 109058, 0, 3,
                                                                       106133, 87548, 106583,
                                                                       66708, 67128, 89888,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 109733, 0, 3,
                                                                       106583, 87908, 107033,
                                                                       67128, 67548, 90428,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 110408, 0, 3,
                                                                       107033, 88268, 107483,
                                                                       67548, 67968, 90968,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 111083, 0, 3,
                                                                       107483, 88628, 107933,
                                                                       67968, 68388, 91508,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 111758, 0, 3,
                                                                       108383, 89348, 109058,
                                                                       69228, 69816, 92048,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 112703, 0, 3,
                                                                       109058, 89888, 109733,
                                                                       69816, 70404, 92804,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 113648, 0, 3,
                                                                       109733, 90428, 110408,
                                                                       70404, 70992, 93560,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 114593, 0, 3,
                                                                       110408, 90968, 111083,
                                                                       70992, 71580, 94316,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 115538, 0, 3,
                                                                       111758, 92048, 112703,
                                                                       72756, 73540, 95072,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 116798, 0, 3,
                                                                       112703, 92804, 113648,
                                                                       73540, 74324, 96080,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 118058, 0, 3,
                                                                       113648, 93560, 114593,
                                                                       74324, 75108, 97088,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 119318, 0, 3,
                                                                       115538, 95072, 116798,
                                                                       76676, 77684, 98096,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 120938, 0, 3,
                                                                       116798, 96080, 118058,
                                                                       77684, 78692, 99392,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 122558, 0, 3,
                                                                       119318, 98096, 120938,
                                                                       80708, 81968, 100688,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 124583, 111758, 945, ncols);

                    simdfunc::contract_primitives(buffer, 125885, 115538, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 127621, 119318, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 129853, 122558, 2025, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 125528, 124583, 21, 1, nmax);

        simdtrf::transform_l_inner(buffer, 127145, 125885, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 129241, 127621, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 131878, 129853, 45, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 132643, 125528, 127145, 17,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 133714, 127145, 129241, 17,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 135142, 129241, 131878, 17,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 136978, 132643, 133714, 17,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 139120, 133714, 135142, 17,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 141976, 136978, 139120, 17,
                                             nmax);

        simdtrf::transform_f_inner(buffer, 145546, 141976, 21, 17, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 145546, 119, nmax);
    }

    for (size_t m = 0; m < 1309; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
